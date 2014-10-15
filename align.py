import numpy as np
import sys

def main():
    (sequence_A, sequence_B, isLocal, gap_penalties, alphabet, match_map) = read_input_file(sys.argv[1])
    #print sequence_A, sequence_B, isLocal, gap_penalties, alphabet, match_map
    alignment(sequence_A, sequence_B, gap_penalties, alphabet, match_map, isLocal)
    #global_alignment(sequence_A, sequence_B, gap_penalties, alphabet, match_map)

    #if(isLocal):
        #do local alignment
    #    print "local alignment"
    #else:
    #    global_alignment(sequence_A, sequence_B, gap_penalties, alphabet, match_map)
        #do global alignment

'''
Take the input file and split it up into its various pieces which are returned to main
'''
def read_input_file(filename):
    content = open(filename).read().splitlines()
    sequence_A = content[0]
    sequence_B = content[1]
    isLocal = False
    if(content[2] == "1"):
        isLocal = True
    gap_penalties = [float(s) for s in content[3].split()]
    alphabet = []
    alphabet.append(content[5])
    alphabet.append(content[7])
    match_map = build_match_map(content[8:])
    return (sequence_A, sequence_B, isLocal, gap_penalties, alphabet, match_map)

'''
Instead of a matrix to lookup characters, I'm building a map
It goes from (letter, letter) -> score
'''
def build_match_map(lines):
    match_map = {}
    for line in lines:
        line_array = [s for s in line.split()]
        if(len(line_array) > 0):
            letter_tuple = (line_array[2], line_array[3])
            match_map[letter_tuple] = float(line_array[4])
    return match_map


def alignment(sequence_A, sequence_B, gap_penalties, alphabet, match_map, isLocal):
    #first we initialize all the variables we need
    width = len(sequence_A) + 1
    height = len(sequence_B) + 1
    M = np.zeros(shape = (height, width))
    Ix = np.zeros(shape = (height, width))
    Iy = np.zeros(shape = (height, width))
    M_pointer_map = {}
    Ix_pointer_map = {}
    Iy_pointer_map = {}

    #then we do the recurrance, as discussed in lecture
    #we start at 1 because row 0 and col 0 are always set = 0
    for i in range(1, height):
        for j in range(1, width):
            letterA = sequence_A[j-1]
            letterB = sequence_B[i-1]

            #The recurrance for M
            M_score = round(M[i-1][j-1] + match_map.get((letterA, letterB)),6)
            Ix_score = round(Ix[i-1][j-1] + match_map.get((letterA, letterB)),6)
            Iy_score = round(Iy[i-1][j-1] + match_map.get((letterA, letterB)),6)
            #we round all the scores to 6 decimal places
            max_score = max(M_score, Ix_score, Iy_score)

            #This is for the backtrace M matrix
            M_pointer_map[i,j] = []
            if(M_score == max_score):
                M_pointer_map[i,j].append("M")
            if(Ix_score == max_score):
                M_pointer_map[i,j].append("Ix")
            if(Iy_score == max_score):
                M_pointer_map[i,j].append("Iy")

            #Recurrance for Ix
            Ix_M_score = M[i-1][j] - gap_penalties[0] #dx (x open)
            Ix_Ix_score = Ix[i-1][j] - gap_penalties[1] #ex (x continue)
            max_Ix_score = max(Ix_M_score, Ix_Ix_score)
            Ix[i][j] = max_Ix_score

            #This is for the backtrace Ix matrix
            Ix_pointer_map[i,j] = []
            if(Ix_M_score == max_Ix_score):
                Ix_pointer_map[i,j].append("M")
            if(Ix_Ix_score == max_Ix_score):
                Ix_pointer_map[i,j].append("Ix")

            #The recurrance for Iy
            Iy_M_score = M[i][j-1] - gap_penalties[2] #dy (y open)
            Iy_Iy_score = Iy[i][j-1] - gap_penalties[3] #ey (y continue)
            max_Iy_score = max(Iy_M_score, Iy_Iy_score)
            Iy[i][j] = max_Iy_score

            #The backtrace for the Iy matrix
            Iy_pointer_map[i,j] = []
            if(Iy_M_score == max_Iy_score):
                Iy_pointer_map[i,j].append("M")
            if(Iy_Iy_score == max_Iy_score):
                Iy_pointer_map[i,j].append("Iy")


            if(isLocal):
                M[i][j] = max(M_score, Ix_score, Iy_score, 0)
            else:
                M[i][j] = max(M_score, Ix_score, Iy_score)

    print M
    print Ix
    print Iy
    print "M traceback"
    print_map_as_matrix(M_pointer_map, width, height)
    print "Ix traceback"
    print_map_as_matrix(Ix_pointer_map, width, height)
    print "Iy traceback"
    print_map_as_matrix(Iy_pointer_map, width, height)
    raw_input("")

    if(isLocal):
        best_score, max_locations = get_local_best_score(M, width, height)
        print "Local: " + str(best_score) + str(max_locations) #WILL I NEED MULTIPLE LOCATIONS??
        alignments = get_different_alignments(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, sequence_A, sequence_B, max_locations, True)
    else:
        best_score, max_locations = get_global_best_score(M, width-1, height-1)
        print "Global: " + str(best_score) + str(max_locations) #This -1,-1 thing is strange
        alignments = get_different_alignments(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, sequence_A, sequence_B, max_locations, False)


#return the alignments as a list of tuples [(ATGC, AGGC), (AATG_C, A__GGC), ....]
#For local I need M to check if the value == 0
def get_different_alignments(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, sequence_A, sequence_B, max_locations, isLocal):
    all_alignments = []

    start = max_locations
    trace_cursor = start

    recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, [], isLocal, trace_cursor, all_alignments)
    print all_alignments
    print len(all_alignments)
    #check if there is a consecutive Ix, Iy or Iy, Ix and delete that.... 
    '''all_alignments = delete_consecutive_Ix_Iy(all_alignments)
    print all_alignments
    print len(all_alignments)
    all_alignments = remove_end_gaps(all_alignments)
    print all_alignments
    print len(all_alignments)
    raw_input("")'''
    readable_alignments = []
    #startB, startA = trace_cursor[0]-1, trace_cursor[1]-1
    #print sequence_A[startA]
    #print sequence_B[startB]
    for i in range(0, len(all_alignments)):
        psuedo_alignment = all_alignments[i]
        #what kind of datastructure to add to 'readable_alignments'
        #(A_sequence, B_sequence)
        A_sequence = ""
        B_sequence = ""
        letter_pointer = trace_cursor[0]-1, trace_cursor[1]-1
        for j in range(0, len(psuedo_alignment)):
            direction = psuedo_alignment[j]
            #if we have M, add (letter, letter)
            if direction == 'M':
                A_sequence += sequence_A[letter_pointer[1]]
                B_sequence += sequence_B[letter_pointer[0]]
                letter_pointer = [letter_pointer[0]-1, letter_pointer[1]-1]
            #if we have Ix, add (letterA, -) (move over one row)
            elif direction == "Ix":
                A_sequence += '_'
                B_sequence += sequence_B[letter_pointer[0]]
                letter_pointer = [letter_pointer[0]-1, letter_pointer[1]]
            #if we have Iy, add (-, letterB)
            elif direction == 'Iy':
                A_sequence += sequence_A[letter_pointer[1]]
                B_sequence += '_'
                letter_pointer = [letter_pointer[0], letter_pointer[1]-1]  
            else:
                print "ERROR 3!!"
        readable_alignments.append((A_sequence[::-1], B_sequence[::-1]))
    print readable_alignments

#removes all the alignments that start or end with Iy or Ix
def remove_end_gaps(all_alignments):
    temp_list = list(all_alignments)
    for alignment in temp_list:
        align_len = len(alignment)
        if(alignment[0] == 'Ix' or alignment[0] == 'Iy'):
            if(alignment in all_alignments):
                all_alignments.remove(alignment)
        if(alignment[align_len-1] == 'Ix' or alignment[align_len-1] == 'Iy'):
            if(alignment in all_alignments):
                all_alignments.remove(alignment)  
    return all_alignments     

#Removes any alignments that have an Ix next to Iy or vice versa
def delete_consecutive_Ix_Iy(all_alignments):
    temp_list = list(all_alignments)
    for alignment in temp_list:
        for i in range(0, len(alignment)-1):
            if(['Ix', 'Iy'] == alignment[i:i+2]):
                if(alignment in all_alignments):
                    all_alignments.remove(alignment)
            elif(['Iy', 'Ix'] == alignment[i:i+2]):
                if(alignment in all_alignments):
                    all_alignments.remove(alignment)
    return all_alignments


#recursively finds all the backtraces
def recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, alignments, isLocal, trace_cursor, all_alignments):
    M_trace = [trace_cursor[0]-1, trace_cursor[1]-1]#Diagonal
    Ix_trace = [trace_cursor[0]-1, trace_cursor[1]] #Going Up (row-1)
    Iy_trace = [trace_cursor[0], trace_cursor[1]-1] #Going Left (col-1)
    if(isLocal):
        if(M[trace_cursor] == 0):
            all_alignments.append(alignments)
            return
    else:
        if(trace_cursor[0] == 0 or trace_cursor[1] == 0): #or == 0
            all_alignments.append(alignments)
            return

    arrows = M_pointer_map.get((trace_cursor[0], trace_cursor[1]))
    alignments.append('M')
    align_M = list(alignments)
    alignments.pop()
    alignments.append('Ix')
    align_Ix = list(alignments)
    alignments.pop()
    alignments.append('Iy')
    align_Iy = list(alignments)
    alignments.pop()

    if(len(arrows) == 3): #always ordered [M, Ix, Iy]
        recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_M, isLocal, M_trace, all_alignments)
        recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Ix, isLocal, Ix_trace, all_alignments)
        recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Iy, isLocal, Iy_trace, all_alignments)
    elif(len(arrows) == 2): #Ordered either [M, Ix], [M, Iy] or [Ix, Iy]
        if(arrows[0] == 'M'):
            recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_M, isLocal, M_trace, all_alignments)
            if(arrows[1] == 'Ix'):
                recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Ix, isLocal, Ix_trace, all_alignments)
            elif(arrows[1] == 'Iy'):
                recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Iy, isLocal, Iy_trace, all_alignments)
        elif(arrows[0] == 'Ix' and arrows[1] == 'Iy'):
            recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Ix, isLocal, Ix_trace, all_alignments)
            recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Iy, isLocal, Iy_trace, all_alignments)
    elif(len(arrows) == 1): #Check for each arrow individually
        if(arrows[0] == 'M'):
            recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_M, isLocal, M_trace, all_alignments)
        elif(arrows[0] == 'Ix'):
            recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Ix, isLocal, Ix_trace, all_alignments)
        elif(arrows[0] == 'Iy'):
            recursive_backtrace(M, M_pointer_map, Ix_pointer_map, Iy_pointer_map, align_Iy, isLocal, Iy_trace, all_alignments)
        else:
            print "ERROR 1!!"
    else:
        print "ERROR 2!!"


#Simple helper function to print maps in a more readable format
def print_map_as_matrix(print_map, width, height):
    for i in range(1, height):
        row = []
        for j in range(1, width):
            row.append(print_map.get((i,j)))
        print row

#Find the best score in the whole M matrix
def get_local_best_score(M, width, height):
    best_score = 0
    best_location = 0
    for i in range(1, height):
        for j in range(1, width):
            temp_score = M[i][j]
            if(temp_score > best_score):
                best_score = temp_score
                best_location = [i, j]
    return (best_score, best_location)

#gets the max_score by going through the last column and last row
def get_global_best_score(M, width, height):
    best_score = 0
    location = [0,0]
    for j in range(1, width):
        temp_score = M[height][j]
        if(temp_score > best_score):
            best_score = temp_score
            location = [height, j]
    for i in range(1, height):
        temp_score = M[i][width]
        if(temp_score > best_score):
            best_score = temp_score
            location = [i, width]
    if(M[height][width] > best_score):
        return(round(M[height, width],1), [height, width])
    return (round(best_score,1), location)




if __name__ == '__main__':
    main()

