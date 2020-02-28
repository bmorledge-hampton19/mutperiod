# This script filters out individual mutations from a bed file as requested through a Tkinter interface.
from TkinterDialog import TkinterDialog, Selections
import os

# This function creates the mutation files with omissions, given the mutation to omit and 
# the original mutation file path.
def omitMutations(mutationFilePath,*mutationsToOmit):
    '''
    This function actually creates the mutation files with omissions, given the mutation to omit and 
    the original mutation file path.
    '''

    print("Preparing to omit", mutationsToOmit, "mutations.")

    mutationsToKeep = list() # This is where we'll store the mutations NOT of the specified omission type.

    # Access the given mutation File
    with open(mutationFilePath, 'r') as mutationFile:

        # Store all of the lines that DON'T have the given mutation
        print("Omitting unwanted mutations.")
        for mutation in mutationFile:
            if mutation.split()[4] not in mutationsToOmit:
                mutationsToKeep.append(mutation)

    # Make the mutation data with omissions directory
    omissionsDirectory = os.path.join(os.path.dirname(mutationFilePath),"mutations_with_omissions")
    if not os.path.exists(omissionsDirectory):
        os.mkdir(omissionsDirectory)

    # Generate a name for the new file with the selected omissions
    mutationsAsText = ""
    for mutation in mutationsToOmit: mutationsAsText += (mutation.replace(">","to")) + "_"
    omissionsFilename = "".join((mutationsAsText,"omitted_",mutationFilePath.rsplit("/",1)[-1]))
    omissionsFilePath = os.path.join(omissionsDirectory,omissionsFilename)

    # Write the stored mutations to the new omissions file.
    print("Writing to new file: ", omissionsFilename)
    with open(omissionsFilePath, 'w') as omissionsFile:
        omissionsFile.writelines(mutationsToKeep)



# Whitespace for AeStHeTiC pUrPoSeS
print()

# Create the list of mutations that can be omitted.
mutations = ("C>A","C>G","C>T","T>A","T>C","T>G")

#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createFileSelector("Bed Mutation File:",0)
#dialog.createDropdown("Mutation to omit:",1,0,options=mutations)
dialog.createLabel("Mutations to Omit:",1,0)
for i,mutation in enumerate(mutations):
    dialog.createCheckbox(mutation, 2+int(i/4), i%4)
dialog.createCheckbox("One file for each omission", 3, 2, 2)
dialog.createReturnButton(4,0,2)
dialog.createQuitButton(4,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
mutationFilePath = list(selections.getFilePaths())[0] # The path to the original bed mutation file
shouldMutationsBeOmitted = list(selections.getToggleStates())[0:6] # A list of the bool values telling what mutations to omit.
createManyFiles = list(selections.getToggleStates())[6] # Should the mutations omitted one at a time, or all together, in one file?
mutationsToOmit = list() # If mutations need to be omitted all at once, we need to keep track of them.

# Send the selected mutations to the omitMutation function to be kicked to the curb.
for i,shouldMutationBeOmitted in enumerate(shouldMutationsBeOmitted):
    if shouldMutationBeOmitted and createManyFiles: omitMutations(mutationFilePath, mutations[i])
    elif shouldMutationBeOmitted: mutationsToOmit.append(mutations[i])

if not createManyFiles: omitMutations(mutationFilePath, *mutationsToOmit)