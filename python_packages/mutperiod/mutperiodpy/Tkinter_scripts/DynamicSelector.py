# A Tk widget wich changes based on the state of a given "commander" widget.

import tkinter as tk
from tkinter import filedialog


class DynamicSelector(tk.Frame):
    """
    A Tk widget which changes based on the state of a given "commander" widget.
    """

    def __init__(self, master, workingDirectory):

        self.workingDirectory = workingDirectory
        self.root = master.root
        super().__init__(master)

        self.controllerSetup = False # Flag to prevent multiple controllers
        self.controllerVar = None # The variable that determines which dynamic display is shown
        self.dynamicDisplays = dict() # The displays corresponding to states of the controller var
        self.currentDisplayKey = None # The current active display.  Helps prevent unnecessary dispaly switching.
        self.defaultDisplayKey = None # Determines which display starts as active.

    def getCurrentDisplay(self): 
        
        # Check to see if we actually have a display for the given key, and initialize an empty display if necessary.
        if not self.currentDisplayKey in self.dynamicDisplays:
            self.initDisplay(self.currentDisplayKey)

        return self.dynamicDisplays[self.currentDisplayKey]

    # Gives me confidence I am accessing Tkinter's strange variable system correctly.
    def getControllerVar(self): return self.controllerVar.get()


    def initDropdownController(self, labelText, options):
        "Initialize a dropdown as the controller for what dynamic content is displayed."

        # Make sure we haven't already initialized another controller
        if self.controllerSetup:
            raise ValueError("Multiple controllers attempted to be set up on one dynamic selector.")
        else: self.controllerSetup = True

        ### Make the dropdown

        # Initialize the stringVar to be modified by the dropdown.
        self.controllerVar = tk.StringVar()
        self.controllerVar.set(options[0])

        # Create a tk.Frame to encompass the dropdown.
        dropdownFrame = tk.Frame(self)
        dropdownFrame.grid(row = 0, column = 0, sticky = tk.W)

        # Add the label
        tk.Label(dropdownFrame, text = labelText).grid(row = 0, column = 0)

        # Initialize the dropdown.
        dropdown = tk.OptionMenu(dropdownFrame, self.controllerVar,*options, command = self.checkController)
        dropdown.grid(row = 0, column = 1, pady = 5, padx = 5, sticky = tk.W)

        # Set the default display variable.
        defaultDisplayKey = options[0]


    def initCheckboxController(self, labelText, default = 0):
        "Initialize a checkbox as the controller for what dynamic content is displayed."

        # Make sure we haven't already initialized another controller
        if self.controllerSetup:
            raise ValueError("Multiple controllers attempted to be set up on one dynamic selector.")
        else: self.controllerSetup = True

        self.controllerVar = tk.IntVar()

        checkbox = tk.Checkbutton(self, text = labelText, variable = self.controllerVar, command = self.checkController)
        checkbox.grid(row = 0, column = 0, pady = 3, sticky = tk.W)

        if default != 0: checkbox.select()


    def initDisplay(self, displayKey, selectionsID = None, row = 1, column = 0, columnSpan = 1, workingDirectory = None):
        "Creates and returns a display corresponding to the given display key, which is a controller variable state."

        from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog

        # Make sure there is no display under this identifier already.
        if displayKey in self.dynamicDisplays:
            raise ValueError("Display key \"" + str(displayKey) + "\" already has an associated display.")

        # Make sure the display does not overlap with the controller var.
        if row == 0 and column == 0: raise ValueError("Display overlaps with controller variable!")

        # Create the dialog object.
        if workingDirectory is None: workingDirectory = self.workingDirectory
        tkinterDialog = TkinterDialog(master = self, title = None, ID = selectionsID, 
                                      row = row, column = column, columnSpan = columnSpan,
                                      workingDirectory = workingDirectory)

        # Add it to the dictionary of dynamic displays and return it.
        self.dynamicDisplays[displayKey] = tkinterDialog
        return tkinterDialog


    def initDisplayState(self):
        "Initializes the dynamic displays by hiding all but the default active display."
        
        # Hide all created displays.
        for displayKey in self.dynamicDisplays:
            self.dynamicDisplays[displayKey].hide()

        # Set the first display
        self.checkController()


    def checkController(self, buffer = None):
        "Called whenever the controller variable is changed to determine what to display"

        # Make sure the dynamic display state has actually changed (or needs initializing.)
        if self.currentDisplayKey is None or self.currentDisplayKey != self.getControllerVar():

            # Hide the old display and enable the currently selected one.
            if self.currentDisplayKey is not None: 
                self.getCurrentDisplay().hide()
            self.currentDisplayKey = self.getControllerVar()
            self.getCurrentDisplay().show()
