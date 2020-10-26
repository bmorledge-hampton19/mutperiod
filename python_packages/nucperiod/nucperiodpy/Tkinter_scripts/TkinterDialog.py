# A modular dialog window used to select relevant files and options for a script.

import tkinter as tk
import tkinter.font as tkFont
import os
from tkinter import filedialog
from typing import List, Dict
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import getFilesInDirectory
from nucperiodpy.Tkinter_scripts.MultipleFileSelector import MultipleFileSelector
from nucperiodpy.Tkinter_scripts.DynamicSelector import DynamicSelector


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


    def createSubDialog(self, row, column = 0, selectionsID = None, columnSpan = 1, workingDirectory = None):
        "Creates and returns another TkinterDialog within the parent dialog."

        # Create the dialog object.
        if workingDirectory is None: workingDirectory = self.workingDirectory
        tkinterDialog = TkinterDialog(master = self, title = None, ID = selectionsID, 
                                      workingDirectory = workingDirectory) 
        tkinterDialog.grid(row = row, column = column, columnspan = columnSpan, sticky = tk.W)

        return tkinterDialog


    def createLabel(self, text, row, column, columnSpan = 1, sticky = True, header = False):
        "Create a simple text label"
        myFont = tkFont.nametofont("TkDefaultFont").copy()
        
        if header: myFont.config(size = "12", weight = "bold")
        label = tk.Label(self, text = text, font = myFont)
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

        # Create a tk.Frame to encompass the dropdown.
        dropdownFrame = tk.Frame(self)
        dropdownFrame.grid(row = row, column = column, columnspan = columnSpan, sticky = tk.W)

        # Add the label
        tk.Label(dropdownFrame, text = labelText).grid(row = 0, column = 0)

        # Initialize the dropdown.
        dropdown = tk.OptionMenu(dropdownFrame, dropdownStringVar,*options)
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
            if dynamicSelector.getCurrentDisplay().ID is not None:
                dynamicSelector.getCurrentDisplay().generateSelections()
                self.selections.addSelections(dynamicSelector.getCurrentDisplay().selections)         

        # Destroy the dialog if this is the top level.
        if self.ID == "Toplevel": self.master.destroy()


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
