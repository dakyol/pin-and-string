import numpy as np

def H(Matrix, R):
    return("")

class FociMatrix():
    def __init__(self):
        self.members = ""

    def create_foci_matrix(self, x_function, y_function, domain):
        matrix =np.array([x_function(domain),y_function(domain)])
        return(matrix)

class FiniteFociMatrix(FociMatrix):
    def __init__(self, x_function, y_function, domain):
        self.matrix = np.array([x_function(domain),y_function(domain)])

    def paired_foci(self, R):
        return(-1)

    def transform(self, R):
        self.transformed_matrix = HTransformationMatrix(self, self.matrix, R)
        return(H(self.matrix,R))

class HTransformationMatrix():
    def __init__(self, FociMatrix, R):

        self.transformed_matrix = H(FociMatrix,R)

class InfiniteFociMatrix(FociMatrix):
    def get_x_function(self):
        return(self.matrix[0])