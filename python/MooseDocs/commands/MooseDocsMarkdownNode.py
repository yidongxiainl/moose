import os
import copy
import bs4
import jinja2
import re
import logging
log = logging.getLogger(__name__)

import MooseDocs
from MooseDocsNode import MooseDocsNode

class MooseDocsMarkdownNode(MooseDocsNode):
  """
  Node for converting
  """
  def __init__(self, md_file=None, syntax=None, parser=None, navigation=None, template=None, template_args=dict(), **kwargs):
    super(MooseDocsMarkdownNode, self).__init__(**kwargs)

    if (not md_file) or (not os.path.exists(md_file)):
      raise Exception('The supplied markdown file must exists: {}'.format(md_file))

    # Extract the MooseLinkDatabase for creating source and doxygen links
    self.__syntax = dict()
    for ext in parser.registeredExtensions:
      if isinstance(ext, MooseDocs.extensions.MooseMarkdown):
        self.__syntax = ext.syntax
        break

    self.__parser = parser
    self.__navigation = navigation
    self.__template = template
    self.__template_args = template_args
    self.__html = None
    self.__md_file = md_file


  def source(self):
    """
    Return the source markdown file.
    """
    return self.__md_file

  def build(self):
    """
    Converts the markdown to html.
    """

    # Read the markdown and parse the HTML
    content, meta = MooseDocs.read_markdown(self.__md_file)
    self.__html = self.__parser.convert(content)

    template_args = copy.copy(self.__template_args)
    template_args['navigation'] = self.__navigation
    template_args.update(meta)

    # Create the template object
    env = jinja2.Environment(loader=jinja2.FileSystemLoader([os.path.join(MooseDocs.MOOSE_DIR, 'docs', 'templates'),
                                                             os.path.join(os.getcwd(), 'templates')]))
    template = env.get_template(self.__template)

    # Render the html via template
    complete = template.render(current=self, **template_args)

    # Make sure the destination directory exists
    destination = self.path()
    if not os.path.exists(destination):
      os.makedirs(destination)

    # Finalize the html
    soup = self.finalize(bs4.BeautifulSoup(complete, 'html.parser'))

    # Write the file
    with open(self.url(), 'w') as fid:
      log.debug('Creating {}'.format(self.url()))
      fid.write(soup.encode('utf-8'))

  def finalize(self, soup):
    """
    Finalize the html:

      1. Converts *.md links to link to the correct html file.
    """

    def finder(node, desired, pages):
      """
      Locate nodes for the 'desired' filename
      """
      if node.source() and node.source().endswith(desired):
        pages.append(node)
      for child in node:
        finder(child, desired, pages)
      return pages

    # Loop over <a> tags and update links containing .md to point to .html
    for link in soup('a'):
      href = link.get('href')
      if href and (not href.startswith('http')) and href.endswith('.md'):
        found = []
        finder(self.root(), href, found)

        # Error if file not found or if multiple files found
        if not found:
          log.error('Failed to locate page for markdown file {} in {}'.format(href, self.source()))
          #link.attrs.pop('href')
          link['class'] = 'moose-bad-link'
          continue

        elif len(found) > 1:
          msg = 'Found multiple pages matching the supplied markdown file {} in {}:'.format(href, self.source())
          for f in found:
            msg += '\n    {}'.format(f.path)
          log.error(msg)

        # Update the link with the located page
        url = self.relpath(found[0].url())
        log.debug('Converting link: {} --> {}'.format(href, url))
        link['href'] = url

    # Fix <pre><code class="python"> to be <pre class="language-python"><code>
    # Add code copy button
    count = 0
    for pre in soup('pre'):
      code = pre.find('code')
      if code and code.has_attr('class'):
        pre['class'] = 'language-{}'.format(code['class'][0])

        if not code.has_attr('id'):
          code['id'] = 'moose-code-block-{}'.format(str(count))
          count += 1

        id = '#{}'.format(code['id'])
        btn = soup.new_tag('button')
        btn['class'] = "moose-copy-button btn"
        btn['data-clipboard-target'] = id
        btn.string = 'copy'

        if "moose-code-div" in pre.parent['class']:
            pre.parent.insert(0, btn)
        else:
            pre.insert(0, btn)

    # Add materialize sections for table-of-contents
    div = soup.find('div', id='moose-markdown-content')
    if div:
      current = div
      for tag in div.contents:
        if isinstance(tag, bs4.element.Tag):
          if tag.name == 'h2':
            current = soup.new_tag('div', id=tag.get('id', '#'))
            current.attrs['class'] = "section scrollspy"

          if current != tag.parent:
            tag.wrap(current)

    # Fix media links
    for img in soup('img'):
      img['src'] = self.relpath(img['src'])

    return soup

  def links(self, repo_url):
    """
    Return the GitHub/GitLab address for editing the markdown file.

    Args:
      repo_url[str]: Web address to use as the base for creating the edit link
    """

    output = [('Edit Markdown', os.path.join(repo_url, 'edit', 'devel', MooseDocs.relpath(self.source())))]

    name = self.breadcrumbs()[-1].name()

    for key, syntax in self.__syntax.iteritems():
      if syntax.hasObject(name):
        include = syntax.filenames(name)[0]
        rel_include = MooseDocs.relpath(include)

        output.append( ('Header', os.path.join(repo_url, 'blob', 'master', rel_include)) )

        source = include.replace('/include/', '/src/').replace('.h', '.C')
        if os.path.exists(source):
          rel_source = MooseDocs.relpath(source)
          output.append( ('Source', os.path.join(repo_url, 'blob', 'master', rel_source)) )

        output.append( ('Doxygen', syntax.doxygen(name)) )

    return output

  def contents(self, level='h2'):
    """
    Return the table of contents.
    """
    soup = bs4.BeautifulSoup(self.content(), 'html.parser')
    for tag in soup.find_all(level):
      if 'id' in tag.attrs:
        yield (tag.contents[0], tag.attrs['id'])

  def content(self):
    """
    Return the generated html from markdown.
    """
    if self.__html is None:
        raise Exception('The "build" command must be executed prior to extracting content.')
    return self.__html

  def url(self, parent=None):
    """
    Return the url to the page to be created.
    """
    path = self.path()
    if parent:
      path = os.path.relpath(parent.path(), path)
    return os.path.join(path, 'index.html')
