# A Widget for a Tkinter dialog that allows for the selection of multiple files.

import os
import tkinter as tk
from tkinter import filedialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getFilesInDirectory


class MultipleFileSelector(tk.Frame):
    "A Widget for a Tkinter dialog that allows for the selection of multiple files."

    def __init__(self, master, title, workingDirectory, fileEnding, additionalFileEndings, *fileTypes):

        # Base class initialization
        super().__init__(master)
        self.grid()

        # Set member variables
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