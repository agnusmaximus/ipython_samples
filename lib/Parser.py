import csv
import scipy.sparse
import numpy as np

class Parser:
    @staticmethod
    def load_matrix(matrix_file, offset=0):
        cr = csv.reader(open(matrix_file))
        rows=[]
        cols=[]
        entries=[]
        for triplet in cr:
            rows.append(int(triplet[0])+offset)
            cols.append(int(triplet[1])+offset)
            entries.append(int(triplet[2]))
        rows = np.array(rows)
        cols = np.array(cols)
        entries = np.array(entries)
        
        return scipy.sparse.csc_matrix((entries,(rows,cols)))
        
    @staticmethod
    def load_dict(dict_file):
        cr = csv.reader(open(dict_file))
        words=[]
        for word_count_pair in cr:
            words.append(word_count_pair[0])
        return words
    @staticmethod
    def load_title(title_file):
        cr = csv.reader(open(title_file))
        titles=[]
        for title_count_pair in cr:
            titles.append(title_count_pair[0])
        return titles
    @staticmethod    
    def load_category(category_file):
        cr = csv.reader(open(category_file))
        categories=[]
        for category_count_pair in cr:
            categories.append(category_count_pair[0])
        return categories

#class Printer:
    