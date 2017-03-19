"""
 File name: word_count.py
 Author: Gordon Kiesling
 VUnetid: Kiesligs
 Email: Gordon.s.kiesling@vanderbilt.edu
 Class: CS2204
 Date: 4/1/16
 Honor statement: I have neither given nor received unauthorized aid on this work
 Assignment Number: 6
 Description: Given two documents of text and a
 list of common words to exclude, find the
 frequency of words in both texts and compare
 the two documents for similarity.
"""

# HELPFUL ADVICE: Attempting to index into a python dictionary using
# a key that is not in the table will throw an error. Before you
# index into a dictionary with a key that did not come from *that*
# dictionary's .keys() method, check that the key is valid.


MINIMUM_APPEARANCE = 0.001
STEM_LENGTH = 6


# main method
def main():
    common_file = "texts/common.txt"
#    fin1 = "texts/atlantis.txt"
#    fin2 = "texts/hamlet.txt"
    fin1 = "texts/sonnets.txt"
    fin2 = "texts/poems.txt"
    fout = "comparison.txt"
    
    # create a set of common words from file
    common_set = read_common_words(common_file)
    # process the two input files, creating a dictionary of work
    atl_dict = process_file(fin1, common_set)
    hamlet_dict = process_file(fin2, common_set)

    print(atl_dict)
    print(hamlet_dict)
    # find the euclidean distance, placing distances in a new dictionary
    # then report the results to the output file
    diff_dic = compare_texts(atl_dict, hamlet_dict)
    report_difference(diff_dic, fout)
    
    # print the vector distance to the console
    # print("Vector Distance:", dist)

    
# read_common_words
# This method reads words from a given file and places them
# into a set, which it returns.
#
# Be aware that in python you could read in some trash from
# each line. You need to sanitize this input (remove non-alpha
# characters. Luckily you have another method to write that
# does just this!
#
# pre:  a valid file name with 1 word per line, words all in lower case
# post: all words in the file are placed in the set & returned
def read_common_words(fname):
    common_set = {}
    file1 = open(fname)
    data = process_word(file1.readline())
    while data != "":
        common_set[data] = 0
        data = process_word(file1.readline())
    file1.close()
    return common_set


# process_file
# This function reads in all words from the given file
# after reading a word it converts it to lower case,
# removes non alphabetic characters and stems it to STEM_LENGTH.
# If the resulting word is considered common or is empty it is ignored.
# Otherwise, the count in the dictionary that matches the word is
# is incremented.
#
# pre:  the name of a text file and a set of words to be ignored.
# post: The file has been read; a dictionary of cleansed words is created and returned
def process_file(fname, common):
    fname_dict = {}
    file = open(fname)
    data = file.read()
    file.close()
    file_list = data.splitlines()

    for x in range(file_list.__len__()):
        word_list = file_list[x].split()
        for y in range(word_list.__len__()):
            word_list[y] = process_word(word_list[y])
            if common.__contains__(word_list[y]) or word_list[y].__len__() == 0:
                pass
            elif fname_dict.__contains__(word_list[y]):
                fname_dict[word_list[y]] += 1
            else:
                fname_dict[word_list[y]] = 1
    return fname_dict


# process_word
# helper function to clean strings up and stem to STEM_LENGTH letters
#
# pre:  a string
# post: the string has been converted to lower case and
#       all non-alphabetic characters have been removed and 
#       the word is reduced to STEM_LENGTH characters.
def process_word(word):
    word = word.lower()
    if not word.isalpha():
        plc = ""
        for x in range(word.__len__()):
            if word[x].isalpha():
                plc += word[x]
        word = plc
    if word.__len__() > STEM_LENGTH:
        word = word[:STEM_LENGTH]
    return word

# compare_texts
# Compares the count dictionaries of 2 texts.
# The result returned is a dictionary of the the square 
# of the euclidean distances for each word with a percentile 
# appearance greater than MINIMUM_APPEARANCE.
#
# pre:  two dictionaries of texts to be compared 
# post: return a new dictionary of words that appeared sufficient
#       times in both input dictionaries. The value associated with
#       each word is its euclidean distance.


def compare_texts(first, second):

    r = max(first.__len__(), second.__len__())
    if first.__len__() == r:
        long = first
        short = second
    else:
        short = first
        long = second

    sum_short = 0
    for key in short:
        sum_short += short[key]
    sum_long = 0
    for key in long:
        sum_long += long[key]
    print(sum_long)
    print(sum_short)

    sam_dic = {}
    for key in long:
        if short.__contains__(key):
            sam_dic[key] = 0
    print(sam_dic)
    print(sam_dic.__len__())

    euc_dict = {}
    for key in long:
        if short.__contains__(key) and long[key]/sum_long > MINIMUM_APPEARANCE \
                and short[key]/sum_short > MINIMUM_APPEARANCE:
            euc_dict[key] = ((long[key]/sum_long) - (short[key]/sum_short))
            euc_dict[key] = euc_dict[key]*euc_dict[key]
    print(euc_dict.__len__())
    return euc_dict


# report_difference
# Accepts a dictionary of shared words and their euclidean distances 
# and computes the total distance between the two texts that produced the dictionary.
#
# This method also prints out the results to a file using the following
# format (start and end tags added for clarity):
# / -- start example output --/
# [ word1, difference1 ]
# [ word2, difference2 ]
# ...
# [ wordN, differenceN ]
# 
# Vector Distance: dist
# /-- end example output --/
#
# word1-wordN are all words whose comparison is include in the final sum.
#   these words are ordered alphabetically.
# difference1-differenceN is the euclidean distance of the word.
# dist is the final total distance, which is also returned.
#
# pre: a dictionary of euclidean distances, and an output filename
# post: the file is filled with the data, and the total distance is returned
def report_difference(diff_dict, fname):
    sorted_keys = sorted(diff_dict)
    output = open(fname, 'w')
    total = 0
    for word in sorted_keys:
        print(diff_dict[word])
        total += diff_dict[word]
        output.write("[ " + word + ", " + repr(round(diff_dict[word], 19)) + " ]" + "\n")
    output.write("\nVector Distance: " + repr(round(total, 14)))
    output.close()
    return total


if __name__ == '__main__':
    main()