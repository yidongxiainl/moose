#!/usr/bin/env python
import os
import traceback
from mooseutils import colorText

try:
    from PyQt5 import QtWidgets, QtCore
    MOOSE_USE_QT5 = True
except:
    MOOSE_USE_QT5 = False

"""
Global for enabling/disabling debug mode.
"""
MOOSE_DEBUG_MODE = False

"""
Global for enabling/disabling testing mode.
"""
MOOSE_TESTING_MODE = False

if MOOSE_USE_QT5:
    class MessageEmitter(QtCore.QObject):
        message = QtCore.pyqtSignal(str, str)

        def write(self, msg, color):
            if not self.signalsBlocked():
                self.message.emit(str(msg), str(color))

    messageEmitter = MessageEmitter()

def mooseMessage(*args, **kwargs):
    """
    A generic message function.

    Args:
        args[tuple]: Comma separated items to be printed, non strings are converted with 'repr'.

    Kwargs:
        error[bool]: (Default: False) When True and 'dialog=True' the the "Critical" icon is included with the message.
        warning[bool]: (Default: False) When True and 'dialog=True' the the "Critical" icon is included with the message.
        traceback[bool]: (Default: False) When True the stack trace is printed with the message.
        dialog[bool]: (Default: False) When true a QDialog object is created, the error will still print to the console as well.
        color[str]: (Default: None) Add the bash color string to the message (see colorText).
        debug[bool]: (Default: False) Print the message only if tools.MOOSE_DEBUG_MODE = True.
        test[bool]: FOR TESTING ONLY! (Default: False) When True the QDialog is not executed, just returned.
        indent[int]: Number of levels to indent (2 spaces are applied to each level)
    """

    # Grab the options
    error = kwargs.pop('error', False)
    warning = kwargs.pop('warning', False)
    trace = kwargs.pop('traceback', False)
    dialog = kwargs.pop('dialog', False)
    color = kwargs.pop('color', None)
    test = kwargs.pop('test', False)
    indent = kwargs.pop('indent', 0)

    # Build the message
    message = []
    for arg in args:
        if not isinstance(arg, str):
            message.append(repr(arg))
        else:
            message.append(arg)
    message = '{}{}'.format(' '*2*indent, ' '.join(message))

    # Show a dialog box
    if MOOSE_USE_QT5 and dialog and not MOOSE_TESTING_MODE:
        box = QtWidgets.QMessageBox()
        box.setText(message)

        if warning:
            box.setIcon(QtWidgets.QMessageBox.Warning)
        elif error:
            box.setIcon(QtWidgets.QMessageBox.Critical)

        if test:
            return box
        box.exec_()

    # Emit the message to any listeners
    messageEmitter.write(message, color)

    # Print the message to screen
    if color:
        message = colorText(message, color)
    print message
    # Show the traceback
    if trace:
        traceback.print_stack()
        stack = ''.join(traceback.format_stack())
        messageEmitter.write(stack, color)


def mooseError(*args, **kwargs):
    """
    A mooseMessage setup to produce an error.
    """
    return mooseMessage('ERROR\n', *args, error = kwargs.pop('error', True),
                                          color = kwargs.pop('color', 'RED'),
                                          traceback = kwargs.pop('traceback', True),
                                          **kwargs)

def mooseWarning(*args, **kwargs):
    """
    A mooseMessage setup to produce a warning.
    """
    return mooseMessage('WARNING\n', *args, warning = kwargs.pop('warning', True),
                                            color = kwargs.pop('color', 'YELLOW'), **kwargs)


def mooseDebug(*args, **kwargs):
    """
    A mooseMessage that only appears with the global MOOSE_DEBUG_MODE = True or debug=True passed directly.

    This automatically indents the message based on how many files exist in the call stack, do remove this simply
    set the 'indent' keyword to zero.
    """
    if kwargs.pop('debug', MOOSE_DEBUG_MODE):
        """
        indent = kwargs.pop('indent', None)
        if not indent:
            files = []
            for item in inspect.stack():
                print item
                if item[1] not in files:
                    files.append(item[1])
            indent = len(files)-3 # -3 removes this file and starts the originating file with 0 indent
        """
        return mooseMessage(*args, color=kwargs.pop('color', 'CYAN'), **kwargs)


def touch(fname):
    """
    Touch a file so to update modified time.
    """
    with open(fname, 'a'):
        os.utime(fname, None)


def gold(filename):
    if not os.path.exists(filename):
        #print mooseError('The supplied filename does not exist:', filename)
        return None

    fn = os.path.basename(filename)
    dn = os.path.dirname(filename)
    gold = os.path.join(dn, 'gold', fn)
    if os.path.exists(gold):
        return gold
    return None

def unique_list(output, input):
    """
    Insert items into list, but only if they are unique.
    """
    for item in input:
        if item not in output:
            output.append(item)
