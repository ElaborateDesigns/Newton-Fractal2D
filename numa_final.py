import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.colors as col
import sys

class fractal2D:
    def __init__(self, pol1, pol2, nrIterations = 10):
        self.pol1 = pol1
        self.pol2 = pol2
        self.nrIterations = nrIterations

    def newtonsMethod(self, guess):
        pol1 = self.pol1
        pol2 = self.pol2
        if pol1 == pol2:
            raise SyntaxError("You have entered the same function")
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
        
        #print(jacobianMatrix)
        for i in range(2): #loop to lambdify the functions in the jacobianMatrix(from sympy to function)
            jacobianMatrix[i,0] = sym.lambdify([a,b], jacobianMatrix[i,0], 'numpy')
            jacobianMatrix[i,1] = sym.lambdify([a,b], jacobianMatrix[i,1], 'numpy')
        print(jacobianMatrix[1,1](1,1))
        
        #start the newton method intervals
        gList = []
        for i in range(self.nrIterations):
            funMatrix[0,0] = pol1(guess[0],guess[1])
            funMatrix[1,0] = pol2(guess[0],guess[1])
            #for i in range(2):
            #    for j in range(2):
            #        jMnumbers[i,j] = jacobianMatrix[i,j](guess[0],guess[1])
            
            jMnumbers[0,0] = jacobianMatrix[0,0](guess[0],guess[1])
            jMnumbers[0,1] = jacobianMatrix[0,1](guess[0],guess[1])
            jMnumbers[1,0] = jacobianMatrix[1,0](guess[0],guess[1])
            jMnumbers[1,1] = jacobianMatrix[1,1](guess[0],guess[1])

            #print(jMnumbers)
            jMnumbers = np.linalg.inv(jMnumbers)
            
            dotProduct = np.dot(jMnumbers,funMatrix)
            guess = guess - dotProduct #xn-1 = xn - J^-1*f(xn)
            gList.append(guess)
            err = 1e-100
            if len(gList) > 2:
                if abs(gList[i][0] - gList[i-1][0]) < err and abs(gList[i][1] - gList[i-1][1]) < err:
                    print(guess)
                    return guess

            print(i)

        print(guess)
        return guess




a = sym.Symbol('a')
b = sym.Symbol('b')
p = fractal2D(a**3 - 3*a*b**2 - 2, 3*a**2*b, 1000)
fractal2D.newtonsMethod(p,[2,3])
