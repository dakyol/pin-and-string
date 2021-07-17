import FociMatrix
import numpy as np

class FiniteFociMatrix(FociMatrix):

    def is_pair(self, R, n1, n2):
        FociMatrix = self.matrix
        [row, column] = FociMatrix.shape

        if column != 2 and row <= 2:
            raise Exception("Matrix must have 2 columns and atleast 3 rows!")
        elif column != 2:
            raise Exception("Matrix must have 2 columns!")
        elif row <= 2:
            raise Exception("Matrix must have atleast 3 rows!")
            

        n1 = np.mod(n1,row)

    #def is_pair(self,R,n1,n2):#self alınca aşaığda kullanacaksak her seferinde matrixi tekrar çağırmaya çalışacağız
        if n1 == n2:
            return False

        leng = self.arc_length(n1,n2)

        M = self.matrix
        [row, column] = M.shape

        x1 = M.take(n1, axis=0, mode='wrap')
        x1_next = M.take(n1+1, axis=0, mode='wrap')
        
        x2 = M.take(n2, axis=0, mode='wrap')
        x2_prev = M.take(n2-1, axis=0, mode='wrap')

        vec_x1 = x1_next - x1
        vec_x2 = x2_prev - x2

        a = vec_x2[0]
        b = -vec_x1[0]

        c = vec_x2[1]
        d = -vec_x1[1]

        M = np.array([[a,b],[c,d]])

        try:
            [alpha, beta] = np.linalg.solve(M,x1-x2)
            y = x1 + alpha*vec_x1
        except np.linalg.LinAlgError as err:
            if 'Singular matrix' in str(err):
                return False #İçten gelenler paralel yapıyor yani olmaz bunlar
            else:
                return("Unknown error occured!")

        if (n1+1-n2)%row == 0:
            cond = True
        else:
            cond = FociMatrix.is_counterclockwise(np.array([x1,y,x2]))

        if cond: #FociMatrix.is_counterclockwise(np.array([x1,y,x2])):
            inner_perimeter = leng + np.abs(alpha)+np.abs(beta)
            inner_condition = R > inner_perimeter
            if inner_condition:
                x1_prev = M.take(n1-1, axis=0, mode='wrap')
                x2_next = M.take(n2+1, axis=0, mode='wrap')
                
                vec_x1 = x1 - x1_prev
                vec_x2 = x2 - x2_next

                a = vec_x2[0]
                b = -vec_x1[0]

                c = vec_x2[1]
                d = -vec_x1[1]

                M = np.array([[a,b],[c,d]])

                try:
                    [alpha, beta] = np.linalg.solve(M,x1-x2)
                    y = x1 + alpha*vec_x1
                except np.linalg.LinAlgError as err:
                    if 'Singular matrix' in str(err):
                        return True #Aynı olmadığını biliyoruz en başta eledik, o zaman paraleller
                    else:
                        return("Unknown error occured!")

                if (n1-1-n2)%row == 0:
                    cond = True
                else:
                    cond = FociMatrix.is_counterclockwise(np.array([x1,y,x2]))

                if cond: #FociMatrix.is_counterclockwise(np.array([x1,y,x2])):
                    if (n1-1-n2)%row == 0:
                        return True #Dışarıdan aldığımız dğrusalsa o zaman bunlar uçtadır yani olur ilki yani içten alınan zaten sağlanyıor
                    else:
                        outer_perimeter = leng + np.abs(alpha)+np.abs(beta)
                        outer_condition = R < outer_perimeter
                        if outer_condition:
                            return True
                        else:
                            return False #Dışarıdan alınan counterclockwise ise R'den küçük olmalı
                else:
                    return True
            else:
                return False #İçeriden alınandan büyük olmalı R
        else:
            return False #they are not pair

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