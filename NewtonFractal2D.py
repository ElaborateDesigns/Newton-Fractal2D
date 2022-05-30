from cmath import nan
from decimal import Overflow
from logging import exception
from math import factorial
from pickle import TRUE
from shutil import register_unpack_format
from turtle import dot
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.colors as col
import sys



###################
### CLASS START ###
###################

    
class fractal2D:
    def __init__(self, pol1, pol2, dxPol1, dyPol1, dxPol2, dyPol2):
        self.pol1 = pol1
        self.pol2 = pol2
        self.dxPol1 = dxPol1
        self.dyPol1 = dyPol1
        self.dxPol2 = dxPol2
        self.dyPol2 = dyPol2

        self.nrIterations = 100
        self.err = 1e-5

        self.zeros = [None]
        self.iterNum = [0,]

        if pol1 == pol2:
            raise SyntaxError("You have entered the same function")
        
        #Grid Calculations stuff
        min_y = 0
        min_x = 0
        max_x = 3
        max_y = 3
        size_x = 4
        size_y = 3
        x = np.linspace (min_x, max_x, size_x)
        y = np.linspace (min_y, max_y, size_y)
        self.grid = np.meshgrid(x,y)
        #print(self.grid)
        self.A = np.empty( (size_x, size_y) )
        self.B = np.zeros(self.A.shape) #default


        
#######################
### CALCULATE ROOTS ###
###     SECTION     ###
#######################        


    def newtonsMethod(self, guess):
        x = guess[0]
        y = guess[1]
        funMatrix = []
        #pol1 = self.pol1
        #pol2 = self.pol2
        #jacobianMatrix = np.array([[a,a],[a,a]]) #the jacobian matrix
        #guess = np.array([[guess[0]],[guess[1]]]) #matrix for the aproximation
        #funMatrix = np.empty([2,1]) #array to input the values given the aproxiamtion in polynom
        #jMnumbers = np.empty([2,2]) #array for the numbers when the jacobian matrix is calculated
        
        
        #jacobianMatrix[0,0] = pol1.diff(a) #dervitive
        #jacobianMatrix[0,1] = pol1.diff(b)
        #jacobianMatrix[1,0] = pol2.diff(a)
        #jacobianMatrix[1,1] = pol2.diff(b)
        
        #pol1 = sym.lambdify([a,b],pol1, 'numpy') #lamdify the polinomes
        #pol2 = sym.lambdify([a,b],pol2, 'numpy')
        #
        #print(jacobianMatrix)
        #for i in range(2): #loop to lambdify the functions in the jacobianMatrix(from sympy to function)
        #    jacobianMatrix[i,0] = sym.lambdify([a,b], jacobianMatrix[i,0], 'numpy')
        #    jacobianMatrix[i,1] = sym.lambdify([a,b], jacobianMatrix[i,1], 'numpy')
        #print(jacobianMatrix[1,1](1,1))
        
        #start the newton method intervals
        for i in range(self.nrIterations):
            funMatrix.append(np.array([self.dxPol1(guess), self.dyPol1(guess)]))
            funMatrix.append(np.array([self.dxPol2(guess), self.dyPol2(guess)]))
            JacobianMatrix.append(funMatrix[0])
            JacobianMatrix.append(funMatrix[1])
            JacobianMatrix = np.array(JacobianMatrix)
        
            #for i in range(2):
            #    for j in range(2):
            #        jMnumbers[i,j] = jacobianMatrix[i,j](guess[0],guess[1])
            

            #print(jMnumbers)
            if np.linalg.det(JacobianMatrix) == 0 and abs(self.pol1(guess)) < 1.e-6 and abs(self.pol2(guess)) < 1.e-9:
                return guess
            elif np.linalg.det(JacobianMatrix) == 0:
                return None
            elif i == len(self.nrIterations)-1:
                return None
            JacobianMatrix=np.array(JacobianMatrix)
            guess = np.array(guess)
            guess = guess - np.dot(np.linalg.inv(JacobianMatrix), np.array([self.pol1([x,y]),self.pol2([x,y])]))
            xNew = guess[0]
            yNew = guess[1]
            if abs(x-xNew) <= 1.e-6 and abs(y-yNew) <= 1.e-6:
                return [xNew, yNew]
            x = xNew
            y = yNew
            #dotProduct = np.dot(jMnumbers,funMatrix)
            #newGuess = guess - dotProduct #xn+1 = xn - J^-1*f(xn)
            ##gList.append(guess)
            #print(guess)
            #print(newGuess)
            #if abs(newGuess[0] - guess[0]) <= 1e-1 and abs(newGuess[1] - guess[1]) <= 1e-1:
            #    print(newGuess)
            #    print(i)
            #    return newGuess
            #if len(gList) > 2:
            #    if abs(gList[i][0] - gList[i-1][0]) < self.err and abs(gList[i][1] - gList[i-1][1]) < self.err:
            #        #print(guess)
            #        return guess

            #print(i)
        

        #print(guess)
        return guess

    
    def simpleNewtonMethod(self, guess):
        self.err = 1e-100
        pol1 = self.pol1
        pol2 = self.pol2
        jacobianMatrix = np.array([[a,a],[a,a]]) #the jacobian matrix
        guess = np.array([[guess[0]],[guess[1]]]) #matrix for the aproximation
        funMatrix = np.empty([2,1]) #array to input the values given the aproxiamtion in polynom
        jMnumbers = np.empty([2,2]) #array for the numbers when the jacobian matrix is calculated
        
        
        
        jacobianMatrix[0,0] = pol1.diff(a) #dervitive
        jacobianMatrix[0,1] = pol1.diff(b)
        jacobianMatrix[1,0] = pol2.diff(a)
        jacobianMatrix[1,1] = pol2.diff(b)
        
        pol1 = sym.lambdify([a,b],pol1, 'numpy') #lamdify the polinomes
        pol2 = sym.lambdify([a,b],pol2, 'numpy')
        

        jacobianMatrix[0,0] = sym.lambdify([a,b], jacobianMatrix[0,0], 'numpy')
        jacobianMatrix[0,1] = sym.lambdify([a,b], jacobianMatrix[0,1], 'numpy')
        jacobianMatrix[1,0] = sym.lambdify([a,b], jacobianMatrix[1,0], 'numpy')
        jacobianMatrix[1,1] = sym.lambdify([a,b], jacobianMatrix[1,1], 'numpy')
        #print(jacobianMatrix)
        #print(jacobianMatrix[1,1](1,1))
        
        #start the newton method intervals
        gList = []
        jMnumbers[0,0] = jacobianMatrix[0,0](guess[0],guess[1])
        jMnumbers[0,1] = jacobianMatrix[0,1](guess[0],guess[1])
        jMnumbers[1,0] = jacobianMatrix[1,0](guess[0],guess[1])
        jMnumbers[1,1] = jacobianMatrix[1,1](guess[0],guess[1])
        #print(guess)
        print(jMnumbers)
        if np.linalg.det(jMnumbers) == 0:
            jMnumbers[0,0] = jMnumbers[0,0] + self.err
            jMnumbers[1,1] = jMnumbers[0,1] + self.err
        jMnumbers = np.linalg.inv(jMnumbers)
        #print(guess)
        for i in range(self.nrIterations):
            funMatrix[0,0] = pol1(guess[0],guess[1])
            funMatrix[1,0] = pol2(guess[0],guess[1])
            #print(guess)
            #for i in range(2):
            #    for j in range(2):
            #        jMnumbers[i,j] = jacobianMatrix[i,j](guess[0],guess[1])
            
            
            dotProduct = np.dot(jMnumbers,funMatrix)
            newGuess = guess - dotProduct #xn-1 = xn - J^-1*f(xn)
            gList.append(guess)
            #try:
            #    guess
            #except RuntimeWarning:
            #    print("something")
            if abs(newGuess[0] - guess[0]) <= 1e-6 and abs(newGuess[1] - guess[1]) <= 1e-6:
                return newGuess
            #if len(gList) > 2:
            #    if abs(gList[i][0] - gList[i-1][0]) < self.err and abs(gList[i][1] - gList[i-1][1]) < self.err:
            #        print(guess)
            #        #print(i)
            #        return guess
        if np.isnan(guess).any() == True:
            print("diverged")
            return nan
        else:
            print(guess[0][0])
            return guess
            
    def getAmountZeros(self,x0):
            
        Zeros=[]
        Tol = 1e-16
        print(len(x0))
        if len(x0) > 1:
            for i in range(len(x0)):
                x1 = [x0[i][0], x0[i][1]]
                #print(i)
                #print(x1)
                q = self.newtonsMethod(x1)
                
                if i == 0: 
                    Zeros.append(q)
                if len(Zeros) > 0:
                    #newZero = False
                #    if not q in Zeros:
                #        Zeros.append(q)
                #else:
                #    Zeros.append(q)

                    for c in range(len(Zeros)):
                        newZero = False
                        #if len(Zeros) > 0:
                        #    if not q in Zeros:
                        #        Zeros.append(q)
                        #else:
                        #    print("ok")
                        
                        t,u = abs(q[0] - Zeros[c][0]),abs(q[1] - Zeros[c][1])
                        #print(t,u)
                        #print(f"This is {t}")
                        #if np.isclose(Zeros[c][0], q[0], atol = Tol) and np.isclose(Zeros[c][1], q[1], atol = Tol):
                        l = np.isclose(Zeros[c][0], q[0], atol = Tol)
                        m = np.isclose(Zeros[c][1], q[1], atol = Tol)
                        if l[0] == True and m[0] == True:
                            newZero = False
                            break
                        else:
                            #print("true")
                            newZero = True
                        #try:
                        #    np.where(np.isclose(Zeros, q, atol = Tol))
                        #if t < Tol and u < Tol:
                        #    #Zeros[c][0] = q[0]
                        #    #Zeros[c][1] = q[1]
                        #    newZero = False
                        #    nzList.append(newZero)
                        #    break
                        #    #print(t,u)
                        #else:
                        #    #Zeros.append(q)
                        #    
                        #    #print(c)
                        #    newZero = True
                        #    nzList.append(newZero)
                        #    #absZero = True
                        #    #break
                        #    #print("same")
                            

                    if newZero == True: #and absZero == False:
                    #if len(nzList) == len(Zeros):
                        Zeros.append(q)
                        print("start")
                        print(len(Zeros))
                        print(i)
                        print(q)
                        print(c)
                        print("end")
                        #break
                #else:
                #    Zeros.append(q)
        else:
            q = self.newtonsMethod(x0)
            Zeros.append(q)
        #for c in range(len(Zeros)):
        #    if c > 1:
        #        try np.where(np.isclose(Zeros, ))
        #        t,u = abs(Zeros[c-1][0] - Zeros[c][0]),abs(Zeros[c-1][1] - Zeros[c][1])
        #        print(f"this is t,u: {t,u}")
        #        if t < 1. and u < 1.:
        #            Zeros[c-1][0] = Zeros[c][0]
        #            Zeros[c-1][1] = Zeros[c][1]
        print(Zeros)
        print(len(Zeros))
        self.Zeros = Zeros
        #print(Zeros)
        #return len(Zeros)
        
        
    
    #def plot(self):
    #    N = 100
    #    xmin, ymin = -10, -10
    #    xmax, ymax = 10, 10
    #    nx, ny = (4, 3)
    #    x = np.linspace(xmin, xmax, nx)
    #    y = np.linspace(ymin, ymax, ny)
    #    xv, yv = np.meshgrid(x, y)
    #    self.getAmountZeros(yv)
        


    
    #def plot2(self): #This is bad and inefficient but I just wanted a working prototype.
#
    #    xmin, ymin = -10, -10
    #    xmax, ymax = 10, 10
    #    nx, ny = (10, 10)
    #    x = np.linspace(xmin, xmax, nx)
    #    y = np.linspace(ymin, ymax, ny)
    #    xv, yv = np.meshgrid(x, y)
    #    #self.getAmountZeros(yv)
    #    self.Zeros=[]
    #    self.A = np.zeros( (nx, ny), dtype=int )
    #    self.B = np.zeros(self.A.shape)
    #    
    #    print(f"Matrix A: \n{self.A}")
    #    x_coord = -1
    #    for x0 in x:
    #        x_coord += 1
    #        y_coord = -1
    #        for y0 in y:
    #            y_coord += 1
    #            #print(f"{x_coord},{y_coord} = {x0},{y0}")
    #            result = self.newtonsMethod( (x0,y0) )
    #            i_coord = -1
    #            found=False
    #            for i in self.Zeros:
    #                i_coord += 1
    #                if (round(float(i[0]),5) == round(float(result[0]),5) 
    #                    and round(float(i[1]),5) == round(float(result[1]),5)  ):
    #                    found=True
    #                    self.A[x_coord,y_coord] = i_coord                         
    #                    break
    #            if (not found): 
    #                self.Zeros.append(result)
    #                self.A[x_coord,y_coord] = i_coord + 1
    #                print( f"New root: {result} ({round(float(result[0]),5)},{round(float(result[1]),5)})")
    #            
    #        
    #    print(f"MATRIX A: \n{self.A}")
    #    self.graph()
                            
                        

    #    self.getAmountZeros(A)
        #self.getAmountZeros(yv)
    def plot2(self, simple = True): #This is bad and inefficient but I just wanted a working prototype. 
        #method = self.newtonsMethod()
        if simple == True:
            method = self.simpleNewtonMethod
        else:
            method = self.newtonsMethod
        
        xmin, ymin = -2, -2
        xmax, ymax = 2, 2
        nx, ny = (10, 10)
        
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        #xv, yv = np.meshgrid(x, y)
        ##self.getAmountZeros(yv)
        self.Zeros=[]
        self.A = np.zeros( (nx, ny), dtype=int )
        self.B = np.zeros(self.A.shape)
        
        print(f"Matrix A: \n{self.A}")
        x_coord = -1
        for x0 in x:
            x_coord += 1
            y_coord = -1
            for y0 in y:
                y_coord += 1
                #print(f"{x_coord},{y_coord} = {x0},{y0}")
                result = method( (x0,y0) )
                i_coord = -1
                found=False
                for i in self.Zeros:
                    i_coord += 1
                    if (round(float(i[0]),5) == round(float(result[0]),5) 
                        and round(float(i[1]),5) == round(float(result[1]),5)  ):
                        found=True
                        self.A[x_coord,y_coord] = i_coord                         
                        break
                if (not found): 
                    self.Zeros.append(result)
                    self.A[x_coord,y_coord] = i_coord + 1
                    print( f"New root: {result} ({round(float(result[0]),5)},{round(float(result[1]),5)})")
                
            
        print(f"MATRIX A: \n{self.A}")
        self.graph()
        
#####################
### PRINT-A-GRAPH ###
###    SECTION    ###
#####################


    def shade_color(self, base, val, minimum, maximum):
      """Shades a base color darker, the more iterations the color is assigned to show."""
      #print(f"shade_color(val={val}, min={min} ")
      #print(f"divider = ( {val}-{min}+1 ")
      if (val < minimum): raise Exception(f"val={val} < minimum={minimum}")
      return tuple(i/(val - minimum + 1) for i in base)

  	
    def make_cmap(self): 
      """Constructs a custom color map according to the number of base colors and shades of these base colors needed"""
      nr_of_colors = 1 + int(np.max(self.A))
      max_iter = int( np.max(self.B) )
      min_iter = int( np.min(self.B) )
      palette = [(2,0.5,0.5), (0.5,2.0,0.5), (0.5,0.5,2), (3,3,0), (0, 2, 2), (2, 0, 2), (0,0,0), (1,1,1) ]
      colors = []
      print(f"nr of colors: {nr_of_colors}, min_iter: {min_iter}, max_iter: {max_iter}")
      for i in range(0, nr_of_colors): 
        base_color = palette[i]
        for j in range(min_iter, max_iter+1): 
          new_color = self.shade_color(base_color, j, min_iter, max_iter)
          colors.append(new_color)
      print(f"cmap: {colors}")
      return col.LinearSegmentedColormap.from_list("mycmap", colors)

    
    def graph(self):
      """Makes a graph from A, a matrix of root indices, and B, a matrix of iteration counts."""
      max_iter = self.nrIterations
      C = self.A*(max_iter + 1) + self.B
      #print(C)
      plt.imshow(C, cmap=self.make_cmap() )
      plt.show()


#################
### CLASS END ###
#################    

a = sym.Symbol('a')
b = sym.Symbol('b')
p = fractal2D(a**8 - 28*a**6*b**2 + 70*a**4*b**4 + 15*a**4 - 28*a**2*b**6 - 90*a**2*b**2 + b**8 + 15*b**4 - 16, 8*a**7*b -56*a**5*b**3 + 56*a**3*b**5 + 60*a**3*b - 8*a*b**7 - 60*a*b**3, 10000)

###################################################################################Polies#####################################################################################

#a**8 - 28*a**6*b**2 + 70*a**4*b**4 + 15*a**4 - 28*a**2*b**6 - 90*a**2*b**2 + b**8 + 15*b**4 - 16
#8*a**7*b -56*a**5*b**3 + 56*a**3*b**5 + 60*a**3*b - 8*a*b**7 - 60*a*b**3

#a**3 - 3*a*b**2 - 2*a - 2
#3*a*2*b - b**3 - 2*b

# a**3 - 3*a*b**2 -1
# 3*a*2*b - b**3

###############################################################################################################################################################################

#for c in range(100):
#    iList = [c*-1,c+1]
#    tlist.append(iList)
#for c in range(100):
#    iList = [c,c+1]
#    tlist.append(iList)
#for c in range(100):
#    iList = [c*-1,(c+1)*-1]
#    tlist.append(iList)

#print(yv)
#tlist=[]
#      
#xmin, ymin = -200, -200
#xmax, ymax = 200, 200
#nx, ny = (2, 1000)
#x = np.linspace(xmin, xmax, nx)
#y = np.linspace(ymin, ymax, ny)
#xv, yv = np.meshgrid(x, y)
#tlist.append(yv)
#print(yv)
#print(len(tlist[0]))
#fractal2D.simpleNewtonMethod(p, [4,6])
#fractal2D.newtonsMethod(p, [10,20])
#fractal2D.plot2(p)
#fractal2D.getAmountZeros(p,yv)

def F(x):
    return x[0] ** 3 - 3 * x[0] * x[1] ** 2 - 1

def G(x):
    return 3 * x[0] ** 2 * x[1] - x[1] ** 3

def dxF(x):
    return 3 * x[0] ** 2 - 3 * x[1] ** 2

def dyF(x):
    return -6 * x[0] * x[1]

def dxG(x):
    return 6 * x[0] * x[1]

def dyG(x):
    return 3 * x[0] ** 2 - 3 * x[1] ** 2


def H(x):
    return x[0]**3 - 3*x[0]*x[1]**2-2*x[0]-2

def I(x):
    return 3*x[0]**2*x[1]-x[1]**3-2*x[1]

def dxH(x):
    return 3*x[0]**2-3*x[1]**2-2

def dyH(x):
    return -6*x[0]*x[1]

def dxI(x):
    return 6*x[0]*x[1]

def dyI(x):
    return 3*x[0]**2-3*x[1]**2-2

def J(x):
    return x[0]**8 -28*x[0]**6*x[1]**2 +70*x[0]**4*x[1]**4 +15*x[0]**4 -28*x[0]**2*x[1]**6 -90*x[0]**2*x[1]**2 +15*x[1]**4 -16

def K(x):
    return 8*x[0]**7*x[1] -56*x[0]**5*x[1]**3 +56*x[0]**3*x[1]**5 +60*x[0]**3*x[1] -8*x[0]*x[1]**7 -60*x[0]*x[1]**3

def dxJ(x):
    return 8*x[0]**7 +280*x[0]**3*x[1]**4 +60*x[0]**3 -168*x[0]**5*x[1]**2 -56*x[0]*x[1]**6 -180*x[0]*x[1]**2

def dyJ(x):
    return 8*x[1]**7 +60*x[1]**3 +280*x[0]**4*x[1]**3 -168*x[0]**2*x[1]**5 -180*x[0]**2*x[1] -56*x[0]**6*x[1]

def dxK(x):
    return 56*x[1]*x[0]**6 -280*x[1]**3*x[0]**4 +168*x[0]**5*x[0]**2 +180*x[1]*x[0]**2 -8*x[1]**7-60*x[1]**3

def dyK(x):
    return 8*x[0]**7 - 168*x[0]**5*x[1]**2 +280*x[0]**3*x[1]**4 + 60*x[0]**3 -56*x[0]*x[1]**6 -180*x[0]*x[1]**2


p = fractal2D(F, G, dxF, dyF, dxG, dyG)
q = fractal2D(H, I, dxH, dyH, dxI, dyI)
r = fractal2D(J, K, dxJ, dyJ, dxK, dyK)