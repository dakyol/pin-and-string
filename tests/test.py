import numpy as np

def H(Matrix, R):
    return("")

class FociMatrix():
    def __init__(self, x_function, y_function, domain):
        self.matrix = np.array([x_function(domain),y_function(domain)])
    
    def is_clockwise(self):
        summa = 0
        M = self.matrix
        [row, column] = M.shape
        for i in range(row-1):
            summa += (M[0][i+1]-M[0][i])*(M[1][i+1]+M[1][i])
        summa += (M[0][0]-M[0][-1])*(M[1][0]+M[1][-1])
        if summa == 0:
            raise Exception("This closed surface is not two dimensional!")
        elif summa > 0:
            result = True
        else:
            result = False
        return(result)

    def is_counterclockwise(self):
        result = not self.is_clockwise()
        return(result)

    def convex_is_clockwise(self):
        M = self.matrix
        if len(M)<3:
            raise Exception("At least three points needed to be given!")
        else:
            matrix = np.array([[M[0][0],M[0][1],M[0][2]],[M[1][0],M[1][1],M[1][2]]])
        return(self.is_clockwise())

    def make_clockwise(self):
        if self.matrix.is_counterclockwise():
            self.matrix = np.flip(self.matrix)
            
    def make_counterclockwise(self):
        if self.matrix.is_clockwise():
            self.matrix = np.flip(self.matrix)

class FiniteFociMatrix(FociMatrix):

    def paired_foci(self, R):
        [row, column] = self.matrix.shape
        for i in range(row):
            
        M = ""
        N = M
        if N.convex_is_clockwise():
            n = 1
        else:
            n = 1
        return(M)

    def transform(self, R):
        self.transformed_matrix = HTransformedMatrix(self, self.matrix, R)
        return(H(self.matrix,R))

class HTransformedMatrix():
    def __init__(self, FociMatrix, R):

        self.transformed_matrix = H(FociMatrix,R)

class InfiniteFociMatrix(FociMatrix):
    def get_x_function(self):
        return(self.matrix[0])