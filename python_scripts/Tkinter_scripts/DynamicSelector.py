import tkinter as tk
from tkinter import filedialog

class DynamicSelector(tk.Frame):
    "A Tk widget wich changes based on the state of a given \"commander\" widget."

    def __init__(self, master, workingDirectory):
        import Tkinter.TkinterDialog as TkD
        self.toplevel = master.toplevel
        self.master = master
        self.workingDirectory = workingDirectory
        super().__init__(self.master)
        self.grid()

        self.controllerSetup = False # Flag to prevent multiple controllers
        self.controllerVar = None # The variable that determines which dynamic display is shown
        self.dynamicDisplays = dict() # The displays corresponding to states of the controller var
        self.currentDisplayKey = None # The current active display.  Helps prevent unnecessary dispaly switching.
        self.defaultDisplayKey = None # Determines which display starts as active.

    def getCurrentDisplay(self): return self.dynamicDisplays[self.currentDisplayKey]

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

        # Add the label
        tk.Label(self, text = labelText).grid(row = 0, column = 0)

        # Initialize the dropdown.
        dropdown = tk.OptionMenu(self, self.controllerVar,*options, command = self.checkController)
        dropdown.grid(row = 0, column = 1, pady = 5, padx = 5, sticky = tk.W)

        # Set the default display variable.
        defaultDisplayKey = options[0]


    def initCheckboxController(self, labelText):
        "Initialize a checkbox as the controller for what dynamic content is displayed."

        # Make sure we haven't already initialized another controller
        if self.controllerSetup:
            raise ValueError("Multiple controllers attempted to be set up on one dynamic selector.")
        else: self.controllerSetup = True

        self.controllerVar = tk.IntVar()

        checkbox = tk.Checkbutton(self, text = labelText, variable = self.controllerVar, command = self.checkController)
        checkbox.grid(row = 0, column = 0, columnspan = 2, pady = 3, sticky = tk.W)


    def initDisplay(self, displayKey, selectionsID, workingDirectory = None):
        "Creates and returns a display corresponding to the given display key, which is a controller variable state."

        from Tkinter.TkinterDialog import TkinterDialog

        # Make sure there is no display under this identifier already.
        if displayKey in self.dynamicDisplays:
            raise ValueError("Display key \"" + str(displayKey) + "\" already has an associated display.")


        # Create the dialog object.
        if workingDirectory is None: workingDirectory = self.workingDirectory
        tkinterDialog = TkinterDialog(master = self, title = None, ID = selectionsID, 
                                          workingDirectory = workingDirectory) 
        tkinterDialog.grid(row = 1, column = 0, columnspan = 2)

        # Add it to the dictionary of dynamic displays and return it.
        self.dynamicDisplays[displayKey] = tkinterDialog
        return tkinterDialog


    def initDisplayState(self):
        "Initializes the dynamic displays by hiding all but the default active display."
        
        # Hide all created displays.
        for displayKey in self.dynamicDisplays:
            self.dynamicDisplays[displayKey].grid_remove()

        # Set the first display
        self.checkController()


    def checkController(self):
        "Called whenever the controller variable is changed to determine what to display"

        # Make sure the dynamic display state has actually changed (or needs initializing.)
        if self.currentDisplayKey is None or self.currentDisplayKey != self.getControllerVar():

            # Hide the old display and enable the currently selected one.
            if self.currentDisplayKey is not None: 
                self.getCurrentDisplay().grid_remove()
            self.currentDisplayKey = self.getControllerVar()
            self.getCurrentDisplay().grid()
