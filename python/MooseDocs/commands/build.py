import os
import sys
import copy
import multiprocessing
import markdown
import markdown_include
import bs4
import shutil
from distutils.dir_util import copy_tree
import logging
log = logging.getLogger(__name__)

import MooseDocs
from MooseDocs.extensions.MooseMarkdown import MooseMarkdown
from MooseDocsNode import MooseDocsNode
from MooseDocsMarkdownNode import MooseDocsMarkdownNode


def build_options(parser, subparser):
  """
  Command-line options for build command.
  """
  build_parser = subparser.add_parser('build', help='Build the documentation for serving on another system.')
  build_parser.add_argument('--disable-threads', action='store_true', help="Disable threaded building.")
  return build_parser

def make_tree(directory, node, **kwargs):
  """
  Create the tree structure of NavigationNode/MooseDocsMarkdownNode objects
  """
  for p in os.listdir(directory):

    path = os.path.join(directory, p)
    if p in ['index.md', 'index.html']:
      continue

    if os.path.isfile(path) and (path.endswith('.md')):
      name = os.path.basename(path)[:-3]
      child = MooseDocsMarkdownNode(name=name, parent=node, md_file=path, **kwargs)

    elif os.path.isdir(path) and (p not in ['.', '..']):
      name = os.path.basename(path)
      md = os.path.join(path, 'index.md')
      if os.path.exists(md):
        child = MooseDocsMarkdownNode(name=name, parent=node, md_file=md, **kwargs)
      else:
        child = MooseDocsNode(name=name, parent=node, **kwargs)
      make_tree(path, child, **kwargs)

def flat(node):
  """
  Create a flat list of pages for parsing and generation.

  Args:
    node[NavigationNode]: The root node to flatten from
  """
  for child in node:
    if isinstance(child, MooseDocsMarkdownNode):
      yield child
    for c in flat(child):
      yield c


class Builder(object):
  """
  Object for building
  """
  def __init__(self, parser, site_dir, template, template_args, navigation):

    self._site_dir = site_dir

    content_dir = os.path.join(os.getcwd(), 'content')
    kwargs = {'parser': parser,
              'site_dir': self._site_dir,
              'navigation': MooseDocs.yaml_load(navigation),
              'template': template,
              'template_args': template_args}
    self._root = MooseDocsMarkdownNode(name='', md_file=os.path.join(content_dir, 'index.md'), **kwargs)
    make_tree(content_dir, self._root, **kwargs)

    self._pages = [self._root] + list(flat(self._root))

  def __iter__(self):
    """
    Allow direct iteration over pages contained in this object.
    """
    return self._pages.__iter__()

  def build(self, disable_threads=False):
    """
    Build all the pages in parallel.
    """

    if disable_threads:
      for page in self._pages:
        page.build()

    else:
      jobs = []
      for page in self._pages:
        p = multiprocessing.Process(target=page.build)
        p.start()
        jobs.append(p)

      for job in jobs:
        job.join()

      self.copyFiles()

  def copyFiles(self):
    """
    Copy the css/js/fonts/media files for this project.
    """

    def helper(src, dst):
      if not os.path.exists(dst):
        os.makedirs(dst)
      if os.path.exists(src):
        copy_tree(src, dst)

    # Copy js/css/media from MOOSE and current projects
    for from_dir in [os.path.join(MooseDocs.MOOSE_DIR, 'docs'), os.getcwd()]:
      helper(os.path.join(from_dir, 'js'), os.path.join(self._site_dir, 'js'))
      helper(os.path.join(from_dir, 'css'), os.path.join(self._site_dir, 'css'))
      helper(os.path.join(from_dir, 'media'), os.path.join(self._site_dir, 'media'))

def build_site(config_file='moosedocs.yml', disable_threads=False, **kwargs):
  """
  The main build command.
  """

  # Load the YAML configuration file
  config = MooseDocs.load_config(config_file, **kwargs)

  # Create the markdown parser
  extensions, extension_configs = MooseDocs.get_markdown_extensions(config)
  parser = markdown.Markdown(extensions=extensions, extension_configs=extension_configs)

  # Create object for storing pages to be generated
  builder = Builder(parser, config['site_dir'], config['template'], config['template_arguments'], config['navigation'])

  # Create the html
  builder.build(disable_threads=disable_threads)
  return config, parser, builder

def build(*args, **kwargs):
  """
  The main build command.
  """
  build_site(*args, **kwargs)
