# This script filters out individual mutations from a bed file as requested through a Tkinter interface.
from Tkinter.TkinterDialog import TkinterDialog, Selections
from helper_scripts.UsefulFileSystemFunctions import dataDirectory
import os

# This function creates the mutation files with omissions, given the mutation to omit and 
# the original mutation file path.
def filterMutations(mutationFilePath,omit,*mutationsToFilter):
    '''
    This function actually creates the mutation files with omissions, given the mutation to omit and 
    the original mutation file path.
    '''
    if omit: print("Preparing to omit", mutationsToFilter, "mutations.")
    else: print("Preparing to keep", mutationsToFilter, "mutations and omit others.")

    mutationsToKeep = list() # This is where we'll store the mutations NOT of the specified omission type.

    # Access the given mutation File
    with open(mutationFilePath, 'r') as mutationFile:

        # Store all of the lines that satisfy the given conditions
        print("Filtering...")
        for mutation in mutationFile:
            if omit and mutation.split()[4] not in mutationsToFilter:
                mutationsToKeep.append(mutation)
            elif not omit and mutation.split()[4] in mutationsToFilter:
                mutationsToKeep.append(mutation)

    # Make the filtered by mutations directory
    filteredRootDirectory = os.path.join(os.path.dirname(mutationFilePath),"filtered_by_mutations")
    if not os.path.exists(filteredRootDirectory):
        os.mkdir(filteredRootDirectory)

    # Make the subdirectory for this specific filtering.
    mutationsAsText = ""
    for mutation in mutationsToFilter: mutationsAsText += (mutation.replace(">","to")) + "_"

    if omit: filteredSubDirectory = os.path.join(filteredRootDirectory,mutationsAsText+"omitted")
    else: filteredSubDirectory = os.path.join(filteredRootDirectory,"just_"+mutationsAsText[:-1])

    if not os.path.exists(filteredSubDirectory):
        os.mkdir(filteredSubDirectory)

    # Generate a path for the new file with the selected omissions
    if omit: filteredFilename = "".join((mutationsAsText,"omitted_",os.path.split(mutationFilePath)[1]))
    else: filteredFilename = "".join(("just_",mutationsAsText,os.path.split(mutationFilePath)[1]))
    filteredFilePath = os.path.join(filteredSubDirectory,filteredFilename)

    # Write the stored mutations to the new omissions file.
    print("Writing to new file: ", filteredFilename)
    with open(filteredFilePath, 'w') as filteredFile:
        filteredFile.writelines(mutationsToKeep)



# Whitespace for AeStHeTiC pUrPoSeS
print()

# Create the list of mutations that can be omitted.
mutations = ("C>A","C>G","C>T","T>A","T>C","T>G")

#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=dataDirectory)
dialog.createFileSelector("Bed Mutation File:",0,("Bed Files",".bed"))
# DEPRECATED: dialog.createDropdown("Mutation to omit:",1,0,options=mutations)
dialog.createLabel("Mutations:",1,0)
for i,mutation in enumerate(mutations):
    dialog.createCheckbox(mutation, 2+int(i/4), i%4)
dialog.createLabel("Actions:",4,0)
dialog.createCheckbox("Omit selected mutations", 5, 0, 2)
dialog.createCheckbox("Keep selected mutations and omit others", 5, 2, 2)
dialog.createCheckbox("Create one file for each selected mutation", 6, 0, 2)
dialog.createReturnButton(7,0,2)
dialog.createQuitButton(7,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
mutationFilePath = list(selections.getFilePaths())[0] # The path to the original bed mutation file
shouldMutationsBeFiltered = list(selections.getToggleStates())[0:6] # A list of the bool values telling what mutations to filter.
omit = list(selections.getToggleStates())[6] # Should the selected mutations be omitted
keep = list(selections.getToggleStates())[7] # Should the selected mutations be kept, and others omitted.
createManyFiles = list(selections.getToggleStates())[8] # Should the mutations omitted one at a time, or all together, in one file?
mutationsToFilter = list() # If mutations need to be filtered all at once, we need to keep track of them.

if omit == keep: raise ValueError("Error: You must select only one option, omit OR keep.")

print("Working in file",os.path.split(mutationFilePath)[1])

# Send the selected mutations to the filterMutations function to be kicked to the curb.
for i,shouldMutationBeFiltered in enumerate(shouldMutationsBeFiltered):
    if shouldMutationBeFiltered and createManyFiles: filterMutations(mutationFilePath, omit, mutations[i])
    elif shouldMutationBeFiltered: mutationsToFilter.append(mutations[i])

if not createManyFiles: filterMutations(mutationFilePath, omit, *mutationsToFilter)