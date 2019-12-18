#For Flu images primarily
from tkinter.filedialog import askdirectory
import os
import sys
import shutil
import re
import numpy as np
from skimage import io as io
import operator
import glob
import math

note_path_list = []

def get_parent_path(note_path):
        
        parent_path = str.rsplit(note_path, os.sep, 1)[0]

        return parent_path

def load_stack(image_list):

        img_collection = io.ImageCollection(image_list)
        img_stack = io.concatenate_images(img_collection)

        return img_stack

def save_stack(eye_path, img_stack):
         
        stack_names = glob.glob(os.path.join(eye_path,"stack*"))
        message = "\n[+]\tStack Saved:\n\t\t"

        try:
                last_stack_name = str.rsplit(stack_names[-1],os.sep,1)[1] # 1 should be same index as -1. hopefully the except only happens here and not from the next line... but that shouldn't happen as long as the only stack name iis from this program and even the the spplit may work but int would throw a TypeError
                stack_index = int((str.rsplit(last_stack_name, "_")[1]).rstrip(".tif"))
                io.imsave(os.path.join(eye_path,"stack_"+str(stack_index+1)+".tif"), img_stack)
                print(message+os.path.join(eye_path,"stack_"+str(stack_index+1)+".tif"))

        except IndexError:
                io.imsave(os.path.join(eye_path,"stack_0.tif"), img_stack)
                print(message+os.path.join(eye_path,"stack_0.tif"))

def formatter(name):
        nameStrip = name.strip("0123456789-0123456789-0123456789.abcdefghijklmnopqrstuvwxyz")
        nameList = nameStrip.split('_')
        nameHMS = int(nameList[1]) #has hours, min, sec as xx-xx-xx without '-'
        namems = float(nameList[2])*(10**-3)  #milisecond value
        #print(nameHMS, namems)
        nameH = math.floor(nameHMS*(10**-4)) #makes the 6 digit nuber have a decimal after the hours, floors it to the interger to get the hour images were taken
        nameM = math.floor((nameHMS - (nameH*(10**4)))*(10**-2)) # subtracts the hours from the min and sec component of the time then does similar to above to get minutes
        nameS = nameHMS - nameH*(10**4) - nameM*(10**2)
        totalSec = nameH*(60**2) + nameM*60 + nameS + namems #should compute total seconds
        
        if nameH > 12:
                nameHr = nameH - 12
        else:
                nameHr = nameH

        return totalSec, nameHr, nameM, nameS

def compare_notes(new_reference_data, old_reference_data):
        if new_reference_data == old_reference_data:
                return True
        #Field size difference divided by smaller size. Want <10%
        new_field_size = float(new_reference_data[0])
        old_field_size = float(old_reference_data[0])
        field_diff = abs(new_field_size-old_field_size)/min(new_field_size, old_field_size)

        #Slices should be equal
        new_slices = new_reference_data[-1]
        old_slices = old_reference_data[-1]

        #Center coordinates should be within 10 "microns" - whatever um means in note.
        new_x = float(new_reference_data[1])
        new_y = float(new_reference_data[2])
        old_x = float(old_reference_data[1])
        old_y = float(old_reference_data[2])
        center_distance = np.sqrt((new_x-old_x)**2 + (new_y-old_y)**2)

        if field_diff <= 0.1 and new_slices==old_slices and center_distance <= 10:
                return True
        else:
                return False


def sorting(source):
        file_list = os.listdir(source)
        print (source)
        destination_apd1 = os.path.join(source, 'reflectance')
        destination_apd2 = os.path.join(source,'fluorescence')

        for each_file in file_list:

                file_path = os.path.join(source,each_file)
                if os.path.isfile(file_path) is True:
                        if 'APD1' in each_file:
                                if not os.path.exists(destination_apd1):
                                    os.makedirs(destination_apd1)
                                    shutil.move(file_path, destination_apd1)       
                                else:
                                    shutil.move(file_path, destination_apd1)       
                        elif 'APD2' in each_file:
                                if not os.path.exists(destination_apd2):
                                    os.makedirs(destination_apd2)
                                    shutil.move(file_path, destination_apd2)
                                else:
                                    shutil.move(file_path, destination_apd2)
                        elif 'notes.txt' in each_file:
                                if not os.path.exists(os.path.join(get_parent_path(source),"processed_flag.txt")): #check if processing has already occured
                                        note_path_list.append(file_path)
                                else:
                                        continue

                elif os.path.isdir(file_path) is True:
                        if file_path == os.path.join(source,'fluorescence') or file_path == os.path.join(source,'reflectance'):
                                continue
                        elif file_path == os.path.join(source, "OCT"):
                                continue
                        else:
                                new_source = os.path.join(source,each_file)
                                sorting(new_source)

def process_notes():
        
        delimiters = "=", "\n", ","
        regexPattern = '|'.join(map(re.escape,delimiters))

        #not sorting the note_list but if needed:
        #note_path_list.sort(key=operations.itemgetter(1,2)) ---- this will sort on list of lists based on position one and two of the sub listts
        
        old_reference_data = None
        note_path_list_length = len(note_path_list)

        for i, note_path in enumerate(note_path_list):

                image_set_path = get_parent_path(note_path)
                eye_path = get_parent_path(image_set_path)
                fluorescence_path = os.path.join(image_set_path, "fluorescence")
                #try/except in the case where I accidentally started to aquire stack but quit before any images were saved
                try:
                        flo_img_name_list = os.listdir(fluorescence_path)
                except:
                        print("/n/n--------------> ERROR: "+fluorescence_path+" NOT LOADED\n\n")
                        continue
                first_image, last_image = flo_img_name_list[0], flo_img_name_list[-1]
                firstSec, firstHr, firstM, firstS = formatter(first_image)
                lastSec, lastHr, lastM, lastS = formatter(last_image)
                imagingTime = lastSec - firstSec
                
                flo_img_list = [os.path.join(fluorescence_path, img_name) for img_name in flo_img_name_list]

                with open(note_path, 'a+') as note_file:
                        note_file.write("\n\n"+str(imagingTime)+' seconds elapsed')
                        note_file.write('\n\nStarted imaging at: ' + str(firstHr) + ':' + str(firstM) + ':' + str(firstS))
                        note_file.write('\nFinished imaging at: ' + str(lastHr) + ':' + str(lastM) + ':' + str(lastS))
                        note_file.seek(0)
                        note_data = note_file.read()
                print("\nWrote data to:\n\t"+note_path)
                
                split_data = re.split(regexPattern, note_data)
                new_reference_data = split_data[10], split_data[13], split_data[14], split_data[20]
                
                if old_reference_data == None:
                        #make new stack
                        img_stack = load_stack(flo_img_list)

                elif i+1==note_path_list_length:
                        save_stack(eye_path, img_stack) #for the last stack, otherwise it doesn't get saved!
                        with open(os.path.join(eye_path, "processed_flag.txt"), "w") as flag: 
                                flag.write("This folder has been processed, to process again delete this note (and if desired, the stacks made by the stack_maker program")

                elif eye_path != previous_eye_path: #forces saving if the eye_path changes.
                        save_stack(previous_eye_path, img_stack)
                        with open(os.path.join(previous_eye_path, "processed_flag.txt"), "w") as flag: 
                                flag.write("This folder has been processed, to process again delete this note (and if desired, the stacks made by the stack_maker program")

                        img_stack = load_stack(flo_img_list)
                else:
                        same_location_flag = compare_notes(new_reference_data, old_reference_data)
                        if same_location_flag is True:
                                #open stack and concat
                                new_stack = load_stack(flo_img_list)
                                img_stack = np.concatenate((img_stack, new_stack))
                        elif same_location_flag is False:
                                
                                #save current stack with previous loops parent_path
                                save_stack(previous_eye_path, img_stack)
                                #load new stack
                                img_stack = load_stack(flo_img_list)

                #This is so that when a new path is encountered the old img_stakc can be saved. Last in the loop!
                old_reference_data = new_reference_data
                previous_image_set_path = image_set_path
                previous_eye_path = eye_path



def main():
    source = askdirectory() 
    if source == '':
            sys.exit("\n\nExited: No file path selected\n\n")
    sorting(os.path.abspath(source))#the abspath should clean up the mixed file separators
    with open(os.path.join(source, "note_list.txt"), "w") as note_text:  
        for note in note_path_list:
                note_text.write("{}\n".format(note))
    process_notes()
    print('\n\n\n\t\t\t\t---------\n\t\t\t\tCompleted\n\t\t\t\t---------')
main()
