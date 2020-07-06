import tkinter as tk
from tkinter import filedialog
from typing import List, Dict
from UsefulFileSystemFunctions import getFilesInDirectory
import os


class TkinterDialog(tk.Frame):
    "A modular dialog window used to select relevant files and options for a script."

    def __init__(self, master=None, workingDirectory = os.path.dirname(os.path.realpath(__file__)), 
                 title = "tkinter dialog", ID = "Toplevel"):
        
        # Initialize the tkinter dialogue using the parent constructor, 
        # and set references to the working directory and master window
        if master is None: 
            self.master = tk.Tk()
            self.toplevel = self.master
        else: 
            self.master = master
            self.toplevel = self.master.toplevel
        super().__init__(self.master)
        self.workingDirectory = workingDirectory
        self.ID = ID

        # Initialize the grid used to organize UI elements
        self.grid() 
        
        # If a title was provided, it is assumed to be a top-level frame and the title and graphic can be set up.
        if title is not None:
            # Set the title
            self.master.title(title)

            # Put in a nice visual!
            img = tk.PhotoImage(file = os.path.dirname(os.path.realpath(__file__)) + "/test_tube.png")
            self.master.iconphoto(False, img)

        self.entries = list() # A list of entry objects created by the script
        self.toggles = list() # A list of toggle objects created by the script
        self.dropdownVars = list() # A list of stringVars associated with dropdowns
        self.checkboxVars = list() # A list of intVars associated with checkboxes
        self.multipleFileSelectors = list() # A list of MultipleFileSelectors, which contain a list of filePaths.
        self.dynamicSelectors: List[DynamicSelector] = list() # A list of DynamicSelector objects, which contian TkDialogs of their own.
        self.selections: Selections = None # A selections object to be populated at the end of the dialog
        

    def createWidgets(self):
        "Mainly used as an example.  Probably don't run this function; run the functions it calls."
        
        print ("Creating Widgets...")

        # Create two file selectors and store the entry objects that will hold the filenames.
        self.createFileSelector("First File: ", 0)
        self.createFileSelector("Second File: ", 1)

        #Create two buttons, one to print the selected filenames, and one to quit
        self.createReturnButton(2, 0, 2)
        self.createQuitButton(2,2,2)


    def createLabel(self, text, row, column, columnSpan = 1, sticky = True):
        "Create a simple text label"
        label = tk.Label(self, text = text)
        label.grid(row = row, column = column, columnspan = columnSpan)
        if sticky: label.grid(sticky = tk.W)


    def createFileSelector(self, title: str, row: int, *fileTypes, defaultFile = "No file Selected", column = 0,
                           columnSpan = 2, verbose = False, newFile = False, directory = False):
        """
        Creates a file selector for choosing relevant files from the dialog.
        The whole file selector is contained in a single tk.Frame object and has 3 components:
        The title for the selector occupies the first column.
        The selector itself as a tk.Entry object.
        The final column holds the "browse" button.
        """

        if verbose: print ("Creating file selector: " + title)

        fileSelectorFrame = tk.Frame(master = self)
        fileSelectorFrame.grid(row = row, column = column, columnspan = columnSpan, sticky = tk.W)

        #Create the title for the selector.
        tk.Label(fileSelectorFrame,text = title).grid(row = 0, column = 0)

        #Create the selector itself as a tk.entry object.
        textField = tk.Entry(fileSelectorFrame, width = 40)
        textField.grid(row = 0, column = 1, columnspan = 2, pady = 10, padx = 5)
        textField.insert(0, defaultFile)
        textField.xview(len(defaultFile))

        #Create the "browse" button.
        tk.Button(fileSelectorFrame, text = "Browse", command = lambda: self.browseForFile(textField,title,newFile,directory,*fileTypes)).grid(row = 0, column = 3)

        self.entries.append(textField)


    def createMultipleFileSelector(self, title: str, row: int, fileEnding, *fileTypes, additionalFileEndings = list(), columnSpan = 2):
        "Create a file selector which can dynamically add and remove multiple file paths to a given group."

        # Create an instance of the the MultipleFileSelector class, and place it in the dialog at the given row.
        multipleFileSelector = MultipleFileSelector(self, title, self.workingDirectory, 
                                                    fileEnding, additionalFileEndings, *fileTypes)
        multipleFileSelector.grid(row = row, columnspan = columnSpan, sticky = tk.W, pady = 10)

        # Keep track of the file selector so we can access the file paths it contains later.
        self.multipleFileSelectors.append(multipleFileSelector)


    def createDropdown(self, labelText, row, column, options, columnSpan = 1):
        "Create a dropdown menu to select something from a list of options."
        
        # Initialize the stringVar to be modified by the dropdown.
        dropdownStringVar = tk.StringVar()
        dropdownStringVar.set(options[0])
        self.dropdownVars.append(dropdownStringVar)

        # Create a tk.Frame to encompass the checkbox.
        checkboxFrame = tk.Frame(self)
        checkboxFrame.grid(row = row, column = column, columnspan = columnSpan)

        # Add the label
        tk.Label(checkboxFrame, text = labelText).grid(row = 0, column = 0)

        # Initialize the dropdown.
        dropdown = tk.OptionMenu(checkboxFrame, dropdownStringVar,*options)
        dropdown.grid(row = 0, column = 1, pady = 5, padx = 5)

    
    def createCheckbox(self, text, row, column, columnSpan = 1):
        "Create a checkbox that holds a bool input from the user."
        
        # Initialize the intVar to be modified by the checkbox
        checkboxIntVar = tk.IntVar()
        self.checkboxVars.append(checkboxIntVar)

        # Create the checkbox
        tk.Checkbutton(self, text = text, variable = checkboxIntVar).grid(row = row, column = column, columnspan = columnSpan, pady = 3, sticky = tk.W)


    def createButton(self, text, row, column, command, columnSpan = 1):
        "Create a button to execute a given function."
        tk.Button(self, text=text, command=command).grid(row = row, column = column, columnspan = columnSpan)


    def createQuitButton(self, row, column, columnSpan = 1):
        """
        DEPRECATED.  Use create Exit Buttons instead.
        Create a button that exits the python script.
        """
        self.createButton("Quit",row,column,quit, columnSpan=columnSpan)


    def createReturnButton(self, row, column, columnSpan = 1):
        """
        DEPRECATED.  Use create Exit Buttons instead.
        Create a button that returns the user input to the selections object
        """
        self.createButton("Go",row,column,self.generateSelections, columnSpan=columnSpan)


    def createExitButtons(self, row, column, columnSpan = 2):
        "Creates the quit and return buttons in a single tk.Frame that centers them within that frame."
        buttonFrame = tk.Frame(self)
        buttonFrame.grid(row = row, column = column, columnspan = columnSpan, sticky = tk.W+tk.E)
        buttonFrame.grid_columnconfigure((0,1), weight = 1)

        tk.Button(buttonFrame, text = "Go", command = self.generateSelections).grid(row = 0, column = 0)
        tk.Button(buttonFrame, text = "Quit", command = quit).grid(row = 0, column = 1)


    def createTextField(self, labelText, row, column, columnSpan = 2, defaultText = "Type here"):
        "Creates an editable text field."

        # Create an instance of the the text field, and place it in the dialog at the given row.
        textField = tk.Frame(master = self)
        textField.grid()
        textField.grid(row = row, column = column, columnspan = columnSpan, pady = 10, sticky = tk.W)

        # Create the label and text box in the text field object.
        tk.Label(textField, text = labelText).grid(row = 0, column = 0, sticky = tk.W)
        textBox = tk.Entry(textField, width = 20)
        textBox.grid(row = 1, columnspan = 2, pady = 2, padx = 5)
        textBox.insert(0, defaultText)

        self.entries.append(textBox)


    def createDynamicSelector(self, row, column, columnSpan = 1, workingDirectory = None):
        "Creates a dynamic selector object at the given location and returns it so it can be further modified."

        if workingDirectory is None: workingDirectory = self.workingDirectory
        dynamicSelector = DynamicSelector(self, workingDirectory)
        dynamicSelector.grid(row = row, column = column, columnspan = columnSpan, sticky = tk.W)
        self.dynamicSelectors.append(dynamicSelector)
        return dynamicSelector


    def browseForFile(self,textField: tk.Entry, title, newFile, directory, *fileTypes):
        "Opens a UI for selecting a file starting from the working directory."

        fileTypes = fileTypes + (("Any File Type", ".*"),)
        
        if directory:

            filename = filedialog.askdirectory(initialdir = self.workingDirectory, title = title)

        else:
            if not newFile:
                filename = filedialog.askopenfilename(filetypes = fileTypes,
                    initialdir = self.workingDirectory, title = title)
            else:
                filename = filedialog.asksaveasfilename(filetypes = fileTypes,
                    initialdir = self.workingDirectory, title = title)

        if (filename != ""):
            textField.delete(0, tk.END)
            textField.insert(0, filename)
            textField.xview(len(filename))
    

    def generateSelections(self):
        "Populates the selections object with the user's input and then destroys the widget."

        individualFilePaths = list() # A list of the filenames selected with the dialog
        filePathGroups = list() # A list of the groups of filepaths selected through MultipleFileSelectors.
        toggleStates = list() # A list of the states of the toggles in the dialog
        dropdownSelections = list() # A list of the selections from dropdown menus

        # Get all the different Selections-relevant variables from this dialog object
        for entry in self.entries:
            individualFilePaths.append(entry.get())

        for multipleFileSelector in self.multipleFileSelectors:
            filePathGroups.append(multipleFileSelector.getFilePaths())

        for checkboxVar in self.checkboxVars:
            toggleStates.append(checkboxVar.get())

        for stringVar in self.dropdownVars:
            dropdownSelections.append(stringVar.get())

        # Generate the Selections object from the above variables.
        self.selections = Selections(self.ID,individualFilePaths,filePathGroups,toggleStates,dropdownSelections)

        # Add any Selections objects from any dialogs included in dynamic displays
        for dynamicSelector in self.dynamicSelectors:
            dynamicSelector.getCurrentDisplay().generateSelections()
            self.selections.addSelections(dynamicSelector.getCurrentDisplay().selections)         

        # Destroy the dialog if this is the top level.
        if self.ID == "Toplevel": self.master.destroy()


class MultipleFileSelector(tk.Frame):
    "A Widget for a Tkinter dialog that allows for the selection of multiple files."

    def __init__(self, master, title, workingDirectory, fileEnding, additionalFileEndings, *fileTypes):

        # Base class initialization
        super().__init__(master)
        self.toplevel = master.toplevel
        self.grid()

        # Set member variables
        self.master = master
        self.title = title
        self.workingDirectory = workingDirectory
        self.fileTypes = fileTypes
        self.fileEnding = fileEnding # The expected ending for the files the dialog is requesting.  
                                     # (Used when sifting through directories for relevant files)
        self.additionalFileEndings = additionalFileEndings
        self.directories = list()
        self.maxPathWidth = 0

        # Set the indentation for the list of file paths displays.
        self.grid_columnconfigure(0, minsize = 40)
        self.grid_columnconfigure(4, weight = 1)

        self.setupMultipleFileSelector()


    def setupMultipleFileSelector(self):
        "Sets up the title for the file selector and the button used to add files."

        #Create the title for the selector.
        tk.Label(self,text = self.title).grid(row = 0, column = 0, columnspan = 2, sticky = "w")

        #Create the "Add Files" and "Add Directories" buttons.
        tk.Button(self, text = "Add Files", command = lambda: self.browseForPaths(self.title,0,*self.fileTypes)).grid(
            row = 0, column = 2, sticky = "w")
        tk.Button(self, text = "Add Directory", command = lambda: self.browseForPaths(self.title,1)).grid(
            row = 0, column = 3, sticky = "w")


    def browseForPaths(self,title,fileOrDirectory,*fileTypes):
        """
        A modification of the TkinterDialog function, "browseForFile".
        This function allows for the selection of multiple files or a directory
        and adds the result to the list of file path displays.
        fileOrDirectory = 0 for file selection and 1 for directory selection.
        """

        fileTypes = fileTypes + (("Any File Type", ".*"),) # Acceptable file types.

        # Open up a file dialog to select files or a directory.
        if fileOrDirectory == 0:
            paths = filedialog.askopenfilenames(filetypes = fileTypes,
                initialdir = self.workingDirectory, title = title)
        elif fileOrDirectory == 1:
            path = filedialog.askdirectory(initialdir = self.workingDirectory, title = title)

        # Change the working directory to one directory up from wherever the last file was selected (for file selection)
        # or the same directory as where the directory was selected (for directory selection)
        if fileOrDirectory == 0:
            if len(paths) > 0: self.workingDirectory = os.path.join(os.path.dirname(paths[0]),"..")
        elif fileOrDirectory == 1:
            if len(path) > 0: self.workingDirectory = os.path.dirname(path)

        # Add the selected file paths to the list of file path displays (only new paths)
        if fileOrDirectory == 0:
            for path in paths:
                if not path in self.getPaths(): self.addPathDisplay(path)
        elif fileOrDirectory == 1:
            if len(path) > 0 and not path in self.getPaths(): self.addPathDisplay(path)


    def addPathDisplay(self, Path):
        "Adds a PathDisplay object to the UI for the given file path."

        # Create the new object.
        newPathDisplay = PathDisplay(self, Path)
        newPathDisplay.grid(column = 1, columnspan = 4, sticky = tk.W)
        if newPathDisplay.pathLabelWidth > self.maxPathWidth:
            self.maxPathWidth = newPathDisplay.pathLabelWidth
            self.updatePathDisplayLengths(self.maxPathWidth)
        else: newPathDisplay.grid_columnconfigure(0, minsize = self.maxPathWidth + 5)

        # Force an update on the toplevel window.
        # RIP: The last hour and a half of my life spent finding these two lines of code
        self.toplevel.update()
        self.toplevel.geometry("")


    def removePathDisplay(self, pathDisplay):
        "Removes a given path display, and checks if the length of the remaining paths needs to be updated."

        # If this was the largest (or one of the largest) path display objects, find the new largest path display
        # and update all other display lengths accordingly.
        if pathDisplay.pathLabelWidth == self.maxPathWidth:

            pathDisplay.destroy() # Get rid of the old object now so we don't count it when looking for a new max.

            newMax = 0
            for remainingPathDisplay in self.getPaths(returnWholeObjects = True):
                if remainingPathDisplay.pathLabelWidth > newMax: 
                    newMax = remainingPathDisplay.pathLabelWidth

            self.maxPathWidth = newMax
            self.updatePathDisplayLengths(self.maxPathWidth)

        # Otherwise, just remove the given path display
        else: pathDisplay.destroy()
            

    def updatePathDisplayLengths(self, maxPathWidth):
        "Updates the minimum width of the label portion of all the path displays."
        for pathDisplay in self.getPaths(returnWholeObjects = True):
            pathDisplay.grid_columnconfigure(0, minsize = maxPathWidth + 5)


    def getPaths(self, returnWholeObjects = False):
        """
        Retrieves the file paths associated with each of the PathDisplay objects in the file selector.
        Alternatively, the whole objects can be returned instead.
        """
        
        paths = list() # The file paths to be returned.

        # Get all widgets in the file selector.
        children = self.winfo_children()

        # Find the PathDisplay objects in the file selector's widgets, and get their file paths.  Then, return them.
        for child in children:
            if isinstance(child, PathDisplay):
                if returnWholeObjects: paths.append(child)
                else: paths.append(child.path)
        return paths


    def getFilePaths(self):
        """Returns all the file paths associated with the MultipleFileSelector."""
        filePaths = list()

        for path in self.getPaths():
            if os.path.isdir(path):
                filePaths += getFilesInDirectory(path,self.fileEnding,*self.additionalFileEndings)
            else:
                filePaths.append(path)
        
        return list(set(filePaths))


class PathDisplay(tk.Frame):
    "A Tk widget that that displays the end of a file path with a button which deletes itself."

    def __init__(self, master: MultipleFileSelector, path):

        # Base class initialization
        super().__init__(master)
        self.path = path

        # Truncate the file path name to 60 characters, if necessary.
        if len(path) < 60: displayName = path
        else: displayName = "..." + path[-60:]

        # Create the label with the file path name and the button which destroys this object.
        pathLabel = tk.Label(self,text = displayName)
        pathLabel.grid(sticky = tk.W)

        removeButton = tk.Button(self, text = "Remove File", command = lambda: master.removePathDisplay(self))
        removeButton.grid(row = 0,column = 1)

        # Retrieve the path's width.
        pathLabel.update()
        self.pathLabelWidth = pathLabel.winfo_width()


class DynamicSelector(tk.Frame):
    "A Tk widget wich changes based on the state of a given \"commander\" widget."

    def __init__(self, master, workingDirectory):
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

    def getCurrentDisplay(self) -> TkinterDialog: return self.dynamicDisplays[self.currentDisplayKey]

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


class Selections:
    "A data structure to hold the results from the TkinterDialog"

    def __init__(self, ID, individualFilePaths = None, filePathGroups = None, 
                 toggleStates = None, dropdownSelections = None):

        # This is a dictionary for storing a list of input values as lists themselves.  (See key below)
        self.selectionSets: Dict[str, List[List]] = dict()
        self.selectionSets[ID] = list()

        # Populate the selectionSet associated with the given ID.  List indices are given below
        # 0: individual file paths
        # 1: file path groups
        # 2: toggle states
        # 3: dropdown selections
        self.selectionSets[ID].append(individualFilePaths)
        self.selectionSets[ID].append(filePathGroups)
        self.selectionSets[ID].append(toggleStates)
        self.selectionSets[ID].append(dropdownSelections)


    # DEPRECATED: diverts to getIndividualFilePaths
    def getFilePaths(self, ID = "Toplevel") -> list:
        return self.getIndividualFilePaths(ID)

    def getIndividualFilePaths(self, ID = "Toplevel") -> list:
        return self.selectionSets[ID][0]

    def getFilePathGroups(self, ID = "Toplevel") -> list:
        return self.selectionSets[ID][1]

    def getToggleStates(self, ID = "Toplevel") -> list:
        return self.selectionSets[ID][2]

    def getDropdownSelections(self, ID = "Toplevel") -> list:
        return self.selectionSets[ID][3]


    def addSelections(self, newSelections):
        "Combines this Selections object with the given Selections object."

        newSelections: Selections
        for ID in newSelections.selectionSets:
            assert ID not in self.selectionSets, (
                   "ID: " + ID + " is already present in the selection sets.")
            self.selectionSets[ID] = newSelections.selectionSets[ID]
