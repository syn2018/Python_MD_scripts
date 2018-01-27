#
# See RandText.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import random

def read_file_into_word_list(filename):
    """The text file can be opened using the "open" function, which returns
    a file object (named here as "inputFile"):

    inputFile = open(filename, 'r')

    See Lutz & Ascher, Learning Python, section 7.2, for basic information
    on file I/O in Python.

    Use inputFile.read() to extract the full text contained in filename, 
    returned as one long string ("text"), and use text.split() to split the 
    text into a list of individual words ("words").
    Return the list of words from the function.
    """
    pass


def make_prefix_dictionary(filename):
    """make and return the prefix dictionary based on the text in filename,
    where the dictionary maps a pair (2-tuple) of words to a list of words
    which follow that pair in the text.

    (1) read text from filename into words

    (2) create empty prefix dictionary: prefix = {}

    (3a) loop over index i in the range(0, len(words)-2), accessing all
    word triples in the list

    See Lutz & Ascher, Learning Python, section 10.3, for information on
    for loops in Python.

    (3b) if the word pair "blah=(words[i], words[i+1])" is not already in the
    prefix dictionary ("if blah not in prefix"), then create a new empty
    list in the dictionary for that pair ("prefix[blah] = []")

    (3c) whether or not blah was in prefix, append the 3rd word "blee" in 
    the triple to the list associated with the prefix pair 
    ("prefix[blah].append(blee)")
    """
    pass


def make_random_text(prefix, num_words=100):
    """make_random_text(prefix, num_words) generates and returns num_words of
    random text based on the triplets contained in prefix.
    
    (1) choose a random starting pair ("current_pair") using the
    random.choice function from the set of keys in the prefix dictionary
    ("random.choice(prefix.keys())")

    Type "pydoc random.choice" from the command line (or "help(random.choice)"
    from within the interpreter once random has been imported).

    (2) initialize a string ("random_text") by concatenating the two
    words in current_pair with a space in between
    ("current_pair[0] + ' ' + current_pair[1]")
    
    (3a) loop over range(num_words-2) to generate the remaining words in
    the random text

    (3b) check to see if current_pair is not a key in the prefix dictionary
    ("if current_pair not in prefix"), and break if True [since the last two
    words in the input text may not have a suffix]

    See Learning Python section 10.2 for information about "break" in Python.

    (3c) randomly choose a suffix from the list of words associated with
    the current_pair ("random.choice(prefix[current_pair])")

    (3d) concatenate to the existing random_text a space and the newly chosen
    random word

    (3e) set the new current_pair to be a tuple containing the last word
    ("current_pair[1]") of the old current_pair and the newly chosen word

    """
    pass

