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
        self.multipleFileSelectors = list() # A list of MultipleFileSelectors, which contain a list of filePaths.
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
        tk.Button(self, text = "Browse", command = lambda: self.browseForFile(textField,title,*fileTypes)).grid(row = row, column = 3)

        self.entries.append(textField)

    # Create a file selector which can dynamically add and remove multiple file paths to a given group.
    def createMultipleFileSelector(self, title: str, row: int, *fileTypes):
        "Create a file selector which can dynamically add and remove multiple file paths to a given group."

        # Create an instance of the the MultipleFileSelector class, and place it in the dialog at the given row.
        multipleFileSelector = MultipleFileSelector(self, title, self.workingDirectory, *fileTypes)
        multipleFileSelector.grid(row = row, columnspan = 4, sticky = tk.W)

        # Keep track of the file selector so we can access the file paths it contains later.
        self.multipleFileSelectors.append(multipleFileSelector)

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
    def browseForFile(self,textField: tk.Entry, title, *fileTypes):
        "Opens a UI for selecting a file starting from the working directory."

        fileTypes = fileTypes + (("Any File Type", ".*"),)

        filename = filedialog.askopenfilename(filetypes = fileTypes,
            initialdir = self.workingDirectory, title = title)
        if (filename != ""):
            textField.delete(0, tk.END)
            textField.insert(0, filename)
    
    # Populates the selections object with the user's input
    def generateSelections(self):
        "Populates the selections object with the user's input"

        individualFilePaths = list() # A list of the filenames selected with the dialog
        filePathGroups = list() # A list of the groups of filepaths selected through MultipleFileSelectors.
        toggleStates = list() # A list of the states of the toggles in the dialog
        dropdownSelections = list() # A list of the selections from dropdown menus

        for entry in self.entries:
            individualFilePaths.append(entry.get())

        for multipleFileSelector in self.multipleFileSelectors:
            filePathGroups.append(multipleFileSelector.getFilePaths())

        for checkboxVar in self.checkboxVars:
            toggleStates.append(checkboxVar.get())

        for stringVar in self.dropdownVars:
            dropdownSelections.append(stringVar.get())

        self.selections = Selections(individualFilePaths,filePathGroups,toggleStates,dropdownSelections)
        self.master.destroy()


# A Widget for a Tkinter dialog that allows for the selection of multiple files.
class MultipleFileSelector(tk.Frame):
    "A Widget for a Tkinter dialog that allows for the selection of multiple files."

    def __init__(self, master, title, workingDirectory, *fileTypes):

        # Base class initialization
        super().__init__(master)
        self.grid()

        # Set member variables
        self.master = master
        self.title = title
        self.workingDirectory = workingDirectory
        self.fileTypes = fileTypes

        # Set the indentation for the list of file paths displays.
        self.grid_columnconfigure(0, minsize = 40)

        self.setupMultipleFileSelector()

    # Sets up the title for the file selector and the button used to add files.
    def setupMultipleFileSelector(self):
        "Sets up the title for the file selector and the button used to add files."

        #Create the title for the selector.
        tk.Label(self,text = self.title).grid(row = 0, column = 0, columnspan = 2, sticky = "w")

        #Create the "Add Files" button.
        tk.Button(self, text = "Add Files", command = lambda: self.browseForFiles(self.title,*self.fileTypes)).grid(
            row = 0, column = 2, sticky = "w")

    # A modification of the TkinterDialog function, "browseForFile".
    # This function allows for the selection of multiple files and adds new files to the list of file path displays.
    def browseForFiles(self,title,*fileTypes):
        '''
        A modification of the TkinterDialog function, "browseForFile".
        This function allows for the selection of multiple files and adds new files to the list of file path displays.
        '''

        fileTypes = fileTypes + (("Any File Type", ".*"),) # Acceptable file types.

        # Open up a file dialog to select files.
        filePaths = filedialog.askopenfilenames(filetypes = fileTypes,
            initialdir = self.workingDirectory, title = title)

        # Change the working directory to the wherever a file was selected.
        if len(filePaths) > 0: self.workingDirectory = os.path.dirname(filePaths[0])

        # Add the selected file paths to the list of file path displays (only new paths)
        for filePath in filePaths:
            if not filePath in self.getFilePaths(): self.addFilePath(filePath)
    
    # Adds a FilePathDisplay object to the UI for the given file path.
    def addFilePath(self, filePath):
        "Adds a FilePathDisplay object to the UI for the given file path."

        # Create the new object.
        FilePathDisplay(self, filePath).grid(column = 1, columnspan = 3, sticky = tk.W)   

        # Force an update on the toplevel window.
        # RIP: The last hour and a half of my life spent finding these two lines of code
        self.master.master.update()
        self.master.master.geometry("")
    
    # Retrieves the file paths associated with each of the FilePathDisplay objects in the file selector
    def getFilePaths(self):
        "Retrieves the file paths associated with each of the FilePathDisplay objects in the file selector"
        
        filePaths = list() # The file paths to be returned.

        # Get all widgets in the file selector.
        children = self.winfo_children()

        # Find the FilePathDisplay objects in the file selector's widgets, and get their file paths.  Then, return them.
        for child in children:
            if isinstance(child, FilePathDisplay):
                filePaths.append(child.filePath)
        return filePaths


# A Tk widget that that displays the end of a file path with a button which deletes itself.
class FilePathDisplay(tk.Frame):
    "A Tk widget that that displays the end of a file path with a button which deletes itself."

    def __init__(self, master, filePath):

        # Base class initialization
        super().__init__(master)
        self.filePath = filePath

        # Truncate the file path name to 60 characters, if necessary.
        if len(filePath) < 60: displayName = filePath
        else: displayName = "..." + filePath[-60:]

        # Standardize the size of the file path name column, so that the proceeding buttons line up nicely.
        self.grid_columnconfigure(0,minsize = 500)

        # Create the label with the file path name and the button which destroys this object.
        tk.Label(self,text = displayName).grid(sticky = tk.W)
        tk.Button(self, text = "Remove File", command = self.destroy).grid(row = 0,column = 1)


# A data structure to hold the results from the TkinterDialog
class Selections:
    "A data structure to hold the results from the TkinterDialog"

    def __init__(self, individualFilePaths = None, filePathGroups = None, toggleStates = None, dropdownSelections = None):

        self.individualFilePaths = individualFilePaths
        self.filePathGroups = filePathGroups
        self.toggleStates = toggleStates
        self.dropdownSelections = dropdownSelections


    # DEPRECATED: diverts to getIndividualFilePaths
    def getFilePaths(self) -> list:
        return self.getIndividualFilePaths()

    def getIndividualFilePaths(self) -> list:
        return self.individualFilePaths

    def getFilePathGroups(self) -> list:
        return self.filePathGroups

    def getToggleStates(self) -> list:
        return self.toggleStates

    def getDropdownSelections(self) -> list:
        return self.dropdownSelections