import tkinter as tk
from tkinter import filedialog
import os

#A modular dialog window used to select relevant files and options for a script.
class TkinterDialog(tk.Frame):
    "A modular dialog window used to select relevant files for a script."

    def __init__(self, master=None, workingDirectory = os.path.dirname(os.path.realpath(__file__)), title = "tkinter dialog"):
        
        # Initialize the tkinter dialogue using the parent constructor, 
        # and set references to the working directory and master window
        if master is None: self.master = tk.Tk()
        else: self.master = master
        super().__init__(self.master)
        self.workingDirectory = workingDirectory

        # Initialize the grid used to organize UI elements
        self.grid() 
        
        # Set the title
        self.master.title(title)

        # Put in a nice visual!
        img = tk.PhotoImage(file = os.path.dirname(os.path.realpath(__file__)) + "/test_tube.png")
        self.master.iconphoto(False, img)

        self.entries = list() # A list of entry objects created by the script
        self.toggles = list() # A list of toggle objects created by the script
        self.dropdownVars = list() # A list of stringVars associated with dropdowns
        self.checkboxVars = list() # A list of intVars associated with checkboxes
        self.selections = None # A selections object to be populated at the end of the dialog
        
    # Example Function
    def createWidgets(self):
        "Mainly used as an example.  Probably don't run this function; run the functions it calls."
        
        print ("Creating Widgets...")

        # Create two file selectors and store the entry objects that will hold the filenames.
        self.createFileSelector("First File: ", 0)
        self.createFileSelector("Second File: ", 1)

        #Create two buttons, one to print the selected filenames, and one to quit
        self.createReturnButton(2, 0, 2)
        self.createQuitButton(2,2,2)

    # Create a simple text label
    def createLabel(self, text, row, column, columnSpan = 1, sticky = True):
        "Create a simple text label"
        label = tk.Label(self, text = text)
        label.grid(row = row, column = column, columnspan = columnSpan)
        if sticky: label.grid(sticky = tk.W)


    def createFileSelector(self, title: str, row: int, *fileTypes, defaultFile = "No file Selected", verbose = False):
        """
        Creates a file selector for choosing relevant files from the dialog.
        The file selector spans 4 columns:
        The title for the selector occupies the first column.
        The selector itself occupies the second and third columns as a tk.Entry object.
        The final column holds the "browse" button.
        """

        if verbose: print ("Creating file selector: " + title)

        #Create the title for the selector.
        tk.Label(self,text = title).grid(row = row, sticky = tk.W)

        #Create the selector itself as a tk.entry object.
        textField = tk.Entry(self, width = 40)
        textField.grid(row = row, column = 1, columnspan = 2, pady = 10, padx = 5)
        textField.insert(0, defaultFile)

        #Create the "browse" button.
        tk.Button(self, text = "Browse", command = lambda: self.browseForFile(textField,*fileTypes)).grid(row = row, column = 3)

        self.entries.append(textField)

    # Create a dropdown menu to select something from a list of options.
    def createDropdown(self, labelText, row, column, options, columnSpan = 1):
        "Create a dropdown menu to select something from a list of options."
        
        # Initialize the stringVar to be modified by the dropdown.
        dropdownStringVar = tk.StringVar()
        dropdownStringVar.set(options[0])
        self.dropdownVars.append(dropdownStringVar)

        # Add the label
        tk.Label(self, text = labelText).grid(row = row, column = column)

        # Initialize the dropdown.
        dropdown = tk.OptionMenu(self, dropdownStringVar,*options)
        dropdown.grid(row = row, column = column+1, columnspan = columnSpan, pady = 5, padx = 5)

    # Create a checkbox that holds a bool input from the user.
    def createCheckbox(self, text, row, column, columnSpan = 1):
        "Create a checkbox that holds a bool input from the user."
        
        # Initialize the intVar to be modified by the checkbox
        checkboxIntVar = tk.IntVar()
        self.checkboxVars.append(checkboxIntVar)

        # Create the checkbox
        tk.Checkbutton(self, text = text, variable = checkboxIntVar).grid(row = row, column = column, columnspan = columnSpan, pady = 3)


    # Create a button to execute a given function.
    def createButton(self, text, row, column, command, columnSpan = 1):
        "Create a button to execute a given function."
        tk.Button(self, text=text, command=command).grid(row = row, column = column, columnspan = columnSpan)

    # Create a button that exits the python script.
    def createQuitButton(self, row, column, columnSpan = 1):
        "Create a button that exits the python script."
        self.createButton("Quit",row,column,quit, columnSpan=columnSpan)

    # Create a button that returns the user input to the selections object
    def createReturnButton(self, row, column, columnSpan = 1):
        "Create a button that returns the user input to the selections object"
        self.createButton("Go",row,column,self.generateSelections, columnSpan=columnSpan)

    # Opens a UI for selecting a file starting from the working directory.    
    def browseForFile(self,textField: tk.Entry, *fileTypes):
        "Opens a UI for selecting a file starting from the working directory."

        fileTypes = fileTypes + (("Any File Type", ".*"),)

        filename = filedialog.askopenfilename(filetypes = fileTypes,
            initialdir = self.workingDirectory)
        if (filename != ""):
            textField.delete(0, tk.END)
            textField.insert(0, filename)
    
    # Populates the selections object with the user's input
    def generateSelections(self):
        "Populates the selections object with the user's input"

        filePaths = list() # A list of the filenames selected with the dialog
        toggleStates = list() # A list of the states of the toggles in the dialog
        dropdownSelections = list() # A list of the selections from dropdown menus

        for entry in self.entries:
            filePaths.append(entry.get())

        for checkboxVar in self.checkboxVars:
            toggleStates.append(checkboxVar.get())

        for stringVar in self.dropdownVars:
            dropdownSelections.append(stringVar.get())
        
        self.selections = Selections(filePaths,toggleStates,dropdownSelections)
        self.master.destroy()


        

# A data structure to hold the results from the TkinterDialog
class Selections:
    "A data structure to hold the results from the TkinterDialog"

    def __init__(self, filePaths = None, toggleStates = None, dropdownSelections = None):

        self.filePaths = filePaths
        self.toggleStates = toggleStates
        self.dropdownSelections = dropdownSelections

    def getFilePaths(self) -> list:
        return self.filePaths

    def getToggleStates(self) -> list:
        return self.toggleStates

    def getDropdownSelections(self) -> list:
        return self.dropdownSelections