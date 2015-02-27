import scipy.sparse
import math
import Options
import numpy
import Parser

class Iterator:
    def __init__(self, options, M):
        self.M = scipy.sparse.csc_matrix(M)
        self.options = options
    def thresh(self, s_vector, num_entries):
        listOfValues = s_vector.T.todense().tolist()[0]
        indices = self.sort(listOfValues)[1]
        thresholded_vec = numpy.zeros(shape=s_vector.T.shape)
        for i in range(num_entries):
            index = indices[i]
            thresholded_vec[0, index] = listOfValues[index]
        return scipy.sparse.csr_matrix(thresholded_vec).T
            
    def top_k_ind(self, s_vector, num_entries):
        listOfValues = s_vector.T.todense().tolist()[0]
        return self.sort(listOfValues)
        
    def sort(self, vector):
        indices = range(len(vector))
        val_ind = map (lambda x : (abs(vector[x]),x), indices)
        sorted_val_ind = sorted(val_ind)
        sorted_val_ind.reverse()
        sorted_indices = map (lambda x : x[1], sorted_val_ind)
        sorted_values = map (lambda x : x[0], sorted_val_ind)
        return (sorted_values, sorted_indices)
    def single_iteration (self):
        
        M = self.M
        m = self.M.shape[0]
        n = self.M.shape[1]
        k_p = self.options.threshold_m
        k_q = self.options.threshold_n
        tolerance = self.options.tolerance
        max_iteration = self.options.max_iteration
        
        
        ini = self.ini_mat();
        p = ini[0]
        q = ini[1]
        obj0 = float("inf")
        converged = False
        iter = 0
        while not converged:
            p_new = self.thresh((M*q), k_p)
            #print 'value is ' + str(sum(map(lambda x: x*x, p_new.T.todense().tolist()[0])))
            #print 'max p is' + str(max(p_new.todense().tolist()[0]))
            p = p_new
            #print 'value is ' + str(p_new.todense())
            #print 'q thresholding: '
            q_new = self.thresh((p.T * M).T, k_q)
            q = q_new/math.sqrt((q_new.T.dot(q_new)).data[0])
            q_test = (p.T * M)
            #print 'max qt is' + str((p.T * M).T.shape)
            #print 'max q is' + str(max(q_new.todense().tolist()[0]))
            updated_mat = M-p_new*(q_new.T)
            obj1 = updated_mat.data.dot(updated_mat.data)
            if (abs(obj1-obj0) <= tolerance or iter >= max_iteration):
                converged = 1
            iter = iter + 1
            obj0 = obj1
        return (p,q)
    
    def multiple_iterations (self):
        term_inds = []
        doc_inds = []
        term_vals = []
        doc_vals = []
        for i in range(self.options.num_pc):
            print 'handling pc# ' + str(i)
            (p,q) = self.single_iteration()
            (p_val, p_ind) = self.top_k_ind(p, self.options.threshold_m)
            (q_val, q_ind) = self.top_k_ind(q, self.options.threshold_n)
            doc_inds.append(p_ind[0:self.options.threshold_m])
            term_inds.append(q_ind[0:self.options.threshold_n])
            doc_vals.append(p_val[0:self.options.threshold_m])
            term_vals.append(q_val[0:self.options.threshold_n])
            self.remove_cols(q_ind[0:self.options.threshold_n])
            self.remove_rows(p_ind[0:self.options.threshold_m])
            
        return ((doc_vals,doc_inds),(term_vals,term_inds))
    
    @staticmethod
    def run(matrix_path, dict_path, title_path, category_path):
        print 'running'
        matrix = Parser.Parser.load_matrix(matrix_path,-1)
        dict = Parser.Parser.load_dict(dict_path)
        titlelist = Parser.Parser.load_title(title_path)
        categorylist = Parser.Parser.load_title(category_path)
        opt = Options.Options(15, 15, 50, 0.0001, 10)
        it = Iterator(opt, matrix)
        word_pcs = []
        title_pcs = []
        ((p_vals,p_inds),(q_vals,q_inds)) = it.multiple_iterations()
        outputfile = open('./output/output.txt','w')
        for pc, pcv in zip(q_inds, q_vals):
            word_pc = []
            for q_ind, q_val in zip(pc, pcv):
                word_pc.append('%.4f\t%s' % (q_val, dict[q_ind]))
            word_pcs.append(word_pc)
            todel = sorted(pc);
            todel.reverse();
            for q_ind in todel:
                del dict[q_ind]
        for pc, pcv in zip(p_inds, p_vals):
            title_pc = []
            for p_ind,p_val in zip(pc,pcv):
                title_pc.append('%.4f\t%s\t%s' % (p_val,titlelist[p_ind],categorylist[p_ind]))
            title_pcs.append(title_pc)
            todel = sorted(pc);
            todel.reverse();
            for p_ind in todel:
                del titlelist[p_ind]
                del categorylist[p_ind]
        topic_num = 0
        for wpc,tpc in zip(word_pcs,title_pcs):
            outputfile.write( '***** TOPIC %d *****\n' % topic_num)
            for row in wpc:
                outputfile.write(row+'\n')
            outputfile.write('-----\n')
            for row in tpc:
                outputfile.write(row+'\n')
            outputfile.write('\n')
            topic_num = topic_num + 1
        #print word_pcs
        #print title_pcs
        
    def remove_rows(self, inds_to_remove):
        self.M = self.M.T
        self.remove_cols(inds_to_remove)
        self.M = self.M.T
    def remove_cols(self, inds_to_remove):
        inds_to_remove = sorted(inds_to_remove)
        inds_to_remove.reverse()
        for ind in inds_to_remove:
            self.remove_col(ind)
        
    def remove_col(self, ind):
        if ind == 0:
            self.M = self.M[:,ind+1:].tocsc()
        elif ind == (self.M.shape[1]-1):
            self.M = self.M[:,0:ind].tocsc()
        else:
            self.M = scipy.sparse.hstack([self.M[:,0:ind],self.M[:,ind+1:]]).tocsc()
        
    def ini_mat(self):
        p = scipy.sparse.csc_matrix(numpy.ones(shape=(self.M.shape[0],1)))
        q = scipy.sparse.csc_matrix(numpy.ones(shape=(self.M.shape[1],1)))
        p = p/math.sqrt(p.T.dot(p).data[0])
        q = q/math.sqrt(q.T.dot(q).data[0])
        return (p,q)

Iterator.run('./input/sparseMatrix.csv', './input/wordlist.csv', './input/documentTitle.csv', './input/category.csv')