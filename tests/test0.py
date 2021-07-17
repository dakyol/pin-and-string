import numpy as np

class FociMatrix():
    #def __init__(self, x_function, y_function, domain):
    #    self.matrix = np.array([x_function(domain),y_function(domain)])
    
    def __init__(self,Matrix):
        size = Matrix.shape
        if len(size) != 2:
            raise Exception("Given array must be 2 dimesional!")
        else:
            [row, column] = size
        
        if column != 2 and row <= 2:
            raise Exception("Matrix must have 2 columns and atleast 3 rows!")
        elif column != 2:
            raise Exception("Matrix must have 2 columns!")
        elif row <= 2:
            raise Exception("Matrix must have atleast 3 rows!")
        else:
            self.matrix = Matrix

    def is_clockwise(M):
        summa = 0
        [row, column] = M.shape
        for i in range(row-1):
            summa += (M[i+1,0]-M[i,0])*(M[i+1,1]+M[i,1])
        summa += (M[0,0]-M[-1,0])*(M[0,1]+M[-1,1])
        if summa == 0:
            raise Exception("This closed surface is not two dimensional!")
        elif summa > 0:
            result = True
        elif summa < 0:
            result = False
        else:
            raise Exception("Unexpected error occurred!")
        return(result)

    def is_counterclockwise(M):
        result = not FociMatrix.is_clockwise(M)
        return(result)

    def convex_is_clockwise(M):
        if len(M)<3:
            raise Exception("At least three points needed to be given!")
        else:
            matrix = np.array([[M[0][0],M[0][1],M[0][2]],[M[1][0],M[1][1],M[1][2]]])
        return(FociMatrix.is_clockwise(matrix))

    def make_clockwise(M):#Bu aslında çok zor sadece counterclockwise ise çeviriyor şu an Hatalı!
        if FociMatrix.is_counterclockwise(M):
            M = np.flip(M)
            
    def make_counterclockwise(M):#Hatalı!
        if FociMatrix.is_clockwise(M):
            M = np.flip(M)
        
    def arc_length(self,n1,n2):
        M = self.matrix
        indices1 = range(n1,n2)
        indices2 = range(n1+1,n2+1)
        D = M.take(indices1, axis=0, mode='wrap')-M.take(indices2, axis=0, mode='wrap')
        S = np.linalg.norm(np.transpose(D),axis=0)
        leng = np.sum(S)
        return(leng)

    def validator(self):
        return True

    def is_pair(self, R, n1, n2):
        if n1 == n2: 
            return False #Hiçbir nokta kendisiyle eş olamaz!

        FociMatrix = self.matrix
        [row, column] = FociMatrix.shape
        
        n1_prev = np.mod(n1-1,row)
        n1 = np.mod(n1,row)
        n1_next = np.mod(n1+1,row)

        n2_prev = np.mod(n2-1,row)
        n2 = np.mod(n2,row)
        n2_next = np.mod(n2+1,row)

        #Saat yönünün tersine FociMatrix'in belirttiği şekil üzerinde n1 indeksli noktadan n2 indeksli noktaya giden yolun uzunluğu.
        arc_length = self.arc_length(n1,n2)

        #Bu isimlendirmeler bir konvensiyon
        RightPairPrev = FociMatrix[n1_prev,:]
        RightPair = FociMatrix[n1,:] 
        RightPairNext = FociMatrix[n1_next,:]

        LeftPairPrev = FociMatrix[n2_prev,:]
        LeftPair = FociMatrix[n2,:]
        LeftPairNext = FociMatrix[n2_next,:]

        RightToNext = RightPairNext - RightPair #Sağ odaktan bir sonraki noktaya bakan vektör.
        RightToPrev = RightPairPrev - RightPair #Sağ odaktan bir önceki noktaya bakan vektör.

        LeftToNext = LeftPairNext - LeftPair #Sol odaktan bir sonraki noktaya bakan vektör.
        LeftToPrev = LeftPairPrev - LeftPair #Sol odaktan bir önceki noktaya bakan vektör.

        RightToLeft = LeftPair - RightPair #Sağdan sola...

        InnerIntersectionMatrix = np.array([[RightToNext[0],-LeftToPrev[0]], #Aradığımız nota şu, sağ odaktan bir sonraki noktaya giden doğru ve sol odaktan bir önceki noktaya giden doğrunun kesişimi.
                                            [RightToNext[1],-LeftToPrev[1]]])#Yani RightPair + (RightPairNext - RightPair)*alpha = LeftPair + (LeftPairPrev - LeftPair)*beta
                                                                             #Veya RightToNext[0]*alpha - LeftToPrev[0]*beta = RightToLeft[0].

        OuterIntersectionMatrix = np.array([[-RightToPrev[0],LeftToNext[0]], #Aradığımız nota şu, sağ odaktan bir önceki noktaya giden doğru ve sol odaktan bir sonraki noktaya giden doğrunun kesişimi.
                                            [-RightToPrev[1],LeftToPrev[1]]])#Yani RightPair + -(RightPairPrev - RightPair)*alpha = LeftPair + -(LeftPairNext - LeftPair)*beta
                                                                             #Veya -RigtToPrev[0]*alpha + LeftToNext[0]*beta = RightToLeft[0].

        try:
            [alpha, beta] = np.linalg.solve(InnerIntersectionMatrix, RightToLeft)
            InnerIntersectionPoint = RightPair + RightToNext*alpha
            FirstCondition = FociMatrix.is_counterclockwise(np.array([RightPair,InnerIntersectionPoint,LeftPair]))
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                if RightToNext[0]/RightToNext[1] == RightToLeft[0]/RightToLeft[1]: #Çakışık doğrular
                    [alpha, beta] = [0, 0] #Bunu leng + abs alpha + abs beta çevre uzunluğuna eşit olacak şekilde değiştirmemiz gerekiyor
                    FirstCondition = True #Counterclockwise sayabiliriz.
                else: #Paralel ama çakışık olmayan doğrular
                    return False #İçten gelenler paralel olamaz!
            else:
                return("Unknown error occured!")

        if FirstCondition:
            InnerPerimeter = arc_length + abs(alpha) + abs(beta)
            if InnerPerimeter < R:
                try:
                    [alpha, beta] = np.linalg.solve(OuterIntersectionMatrix, RightToLeft)
                    OuterIntersectionPoint = RightPair - RightToPrev*alpha
                    SecondCondition = FociMatrix.is_counterclockwise(np.array([RightPair,OuterIntersectionPoint,LeftPair]))
                except np.linalg.LinAlgError as err:
                    if 'Singular matrix' in str(err):
                        return True #Tersten yapılan paralel veya çakışıksa, R sonlu değilmiş gibi davranırız.
                    else:
                        return("Unknown error occured!")
                if SecondCondition: #Saat yönü kontrolü, kesişim noktası hangi tarafta anlamak için.
                    OuterPerimeter = arc_length + abs(alpha) + abs(beta)
                    if OuterPerimeter > R: #Dışarıdan alınan uzunluk R2den daha kısa olsaydı ilk koşul gibi bu ikisi olmazdı.
                        return True
                    else:
                        return False
                else:
                    return True #ilkini geçtikten sonra buraya geldiyse, arkadan yapılan ters yönde kalsa bile durumu etkilemez.
            else:
                return False
        else:
            return False

    def paired_foci(self, R): #Karmaşıklığı 3n*O("is_pair()") arttırıyor (+3n*k)
        M = self.matrix
        [row, column] = M.shape
        pair_list = -np.ones([row*row,column]) #preallocation

        index = range(row)

        for i in index:
            if self.is_pair(R,0,i):#FiniteFociMatrix.is_pair(M,R,0,i)
                pair_list[0,:] = [0,i]
                righto = i

        #her ikisi de false ise kesin 1 righto+1 eş
        lefto = 0
        index = 1
        while lefto < row-1:#bilmem olabilir
            if self.is_pair(R,lefto+1,righto): #FiniteFociMatrix.is_pair(M,R,lefto+1,righto):
                lefto += 1
            elif self.is_pair(R,lefto,righto+1): #FiniteFociMatrix.is_pair(M,R,lefto,righto+1):
                righto += 1
            else:#FiniteFociMatrix.is_pair(M,R,lefto+1,righto+1)
                lefto += 1
                righto += 1
            pair_list[index,:]=[lefto,righto]
            index += 1
        #index+1:row*row kısmını pair_listin silebilirsin hep -1 olacak zaten sonra returnlersin
        deneme = pair_list
        pair_list = np.delete(pair_list,range(index,row*row),axis=0)
        return([pair_list,deneme])

    def transform(self,R):
        pair_list = self.paired_foci(R)
        M = self.matrix
        x=1