import re
import os
from markdown.util import etree
import logging
log = logging.getLogger(__name__)

import MooseDocs
from markdown.inlinepatterns import Pattern
from MooseCommonExtension import MooseCommonExtension
import utils

class MooseImageFile(MooseCommonExtension, Pattern):
  """
  Markdown extension for handling images.

  Usage:
   !image image_file.png|jpg|etc attribute=setting

  Settings:
    caption[str]: Creates a figcaption tag with the supplied text applied.

  All image filenames should be supplied as relative to the docs directory, i.e., media/my_image.png
  """

  # Find !image /path/to/file attribute=setting
  RE = r'^!image\s+(.*?)(?:$|\s+)(.*)'

  def __init__(self, markdown_instance=None, **kwargs):
    MooseCommonExtension.__init__(self, **kwargs)
    Pattern.__init__(self, self.RE, markdown_instance)
    self._settings['caption'] = None

  def handleMatch(self, match):
    """
    process settings associated with !image markdown
    """

    # A tuple separating specific MOOSE documentation features (self._settings) from HTML styles
    settings = self.getSettings(match.group(3))

    # Read the file and create element
    rel_filename = match.group(2)
    filename = MooseDocs.abspath(match.group(2))
    if not os.path.exists(filename):
      return self.createErrorElement('File not found: {}'.format(rel_filename))

    # Create the figure element
    el = self.applyElementSettings(etree.Element('div'), settings)

    card = etree.SubElement(el, 'div')
    card.set('class', 'card')

    img_card = etree.SubElement(card, 'div')
    img_card.set('class', 'card-image')

    img = etree.SubElement(img_card, 'img')
    img.set('src', os.path.relpath(filename, os.getcwd()))
    img.set('class', 'materialboxed')

    # Add caption
    if settings['caption']:
      caption = etree.SubElement(card, 'div')
      p = etree.SubElement(caption, 'p')
      p.set('class', 'moose-caption')
      p.set('align', "justify")
      p.text = settings['caption']

    return el
