import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Frist we define all functions and their derivatives, dx is derivative with respect to x and dy to y. This includes the one from task 4 and the two from task 8
def F(x): # Since we will send through a list of variables ([x,y]) we can use x[0] as x and x[1] as y
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
    return x[0]**8 -28*x[0]**6*x[1]**2 + 70*x[0]**4*x[1]**4 + 15*x[0]**4 -28*x[0]**2*x[1]**6 - 90*x[0]**2*x[1]**2 + x[1]**8 + 15*x[1]**4 - 16
 
def dxJ(x):
    return 8*x[0]**7 - 6*28*x[0]**5*x[1]**2 + 70*4*x[0]**3*x[1]**4 + 4*15*x[0]**3 - 28*2*x[0]*x[1]**6 - 90*2*x[0]*x[1]**2

def K(x):
    return 8*x[0]**7*x[1] -56*x[0]**5*x[1]**3 +56*x[0]**3*x[1]**5 +60*x[0]**3*x[1] -8*x[0]*x[1]**7 -60*x[0]*x[1]**3
 
def dyJ(x):
    return -28*x[0]**6*2*x[1] + 4*70*x[0]**4*x[1]**3 - 28*6*x[0]**2*x[1]**5 - 2*90*x[0]**2*x[1] + 8*x[1]**7 + 4*15*x[1]**3

def dyJ(x):
    return -28*x[0]**6*2*x[1] + 4*70*x[0]**4*x[1]**3 - 28*6*x[0]**2*x[1]**5 - 2*90*x[0]**2*x[1] + 8*x[1]**7 + 4*15*x[1]**3 

def dxK(x):
    return 7*8*x[0]**6*x[1] - 5*56*x[0]**4*x[1]**3 + 3*56*x[0]**2*x[1]**5 + 3*60*x[0]**2*x[1] - 8*x[1]**7 - 60*x[1]**3

def dyK(x):
    return 8*x[0]**7 - 3*56*x[0]**5*x[1]**2 + 56*5*x[0]**3*x[1]**4 + 60*x[0]**3 - 7*8*x[0]*x[1]**6 - 60*3*x[0]*x[1]**2


class fractal2D:                                                                                    # We create the class fractal2D
    def __init__(self, pol1, pol2, dxPol1 = True, dyPol1 = True, dxPol2 = True, dyPol2 = True):     # Initialized with our polinomial and their derivatives. They are set to True because we want to be able to have them be optional 
        self.pol1 = pol1
        self.pol2 = pol2
        self.dxPol1 = dxPol1
        self.dyPol1 = dyPol1
        self.dxPol2 = dxPol2
        self.dyPol2 = dyPol2

        self.zeroes = [None]            # We make a list called zeroes in which we will store the zeros we find later in the code. The first index "None" indicates divergence
        self.zeroIndex = [0,]           # Here we store the index of iterationsit took to converge, if the value diverges it will just return the first value "0"

        if dxPol1 == True and dyPol1 == True and dxPol2 == True and dyPol2 == True:                 # If no derivatives are given
            self.derivativeAprrox = True                                                            # We set our derivateApprox variable to true which will later use another method to approximate the equivilence of dxpol1 dypol1, and so on 
        else:
            self.derivativeAprrox = False                                                           # If derivatives were given then we do not approximate
    def approxXDer(self, fun, guess):   # Method for derivative approximation in repsect to x
        x = guess[0]
        y = guess[1]
        h = 1.e-6                       # If h is changed so will our approximation
        return (fun([x + h, y]) - fun([x, y])) / h                                                  # Retruns our function with the values x + h and just y since y isnt derivated. (derivatas def)
    def approxYDer(self, fun, guess):   # Method for derivateive with repect to y
        x = guess[0]
        y = guess[1]
        h = 1.e-6
        return (fun([x, y + h]) - fun([x, y])) / h

    def newtonsMethod(self, guess, itersteps = False, combo=False):              # Newton raphson method is given two arguemnts (excluding self), the first is a list, our guess, which is just a point on a 2d plane. And then a boolean iterSteps, this is so that later this can return the amount of iterations it took to converge rather than the zero it converged towards
        x = guess[0]                      # x is given the first value of the guess
        y = guess[1]
        if self.derivativeAprrox == True:                           # If we dont have the functions derivatives then derivativeApprox will be true and our branch of newton will use the approximation instead
            for i in range(100):                                    # The range is the max iteration number so this will continue either until it stops by the returns or go trhough 100 loops 
                JacobianMatrix = [np.array([self.approxXDer(self.pol1, guess), self.approxYDer(self.pol1, guess)]), np.array([self.approxXDer(self.pol2, guess), self.approxYDer(self.pol2, guess)])]   # The jacobian matrix will be the np array of dxpol1, dypol1 and dxpol2, dypol2 in the form of the matrix [[dxpol1, dypol1],[dxpol2, dypol2]]
                JacobianMatrix = np.array(JacobianMatrix)           # The jocaboian is turned into an np.array
                if np.linalg.det(JacobianMatrix) == 0 and abs(self.pol1(guess)) < 1.e-6 and abs(self.pol2(guess)) < 1.e-9:  # If the determinant of the jacobian is 0 we cannot invert it, but if our guess is sufficiently small we can asume that we have reached convergence
                    if combo == True:
                        return (guess, i+1)
                    elif itersteps == False:                          # If itersteps is False we just return the zero we converged towards
                        return guess
                    else:                                           # Else if itersteps if on we are after the amount of iterations it took to converge
                        return i+1
                elif np.linalg.det(JacobianMatrix) == 0:
                    if combo == True:
                        return (None, 0)
                    elif itersteps == False:
                        return None
                    else:
                        return 0
                elif i == 99:                                       # If we run out of range then we conclude that we have diveged
                    if combo == True:
                        return (None, 0)
                    elif itersteps == False:
                        return None
                    else:
                        return 0
                JacobianMatrix = np.array(JacobianMatrix)
                guess = np.array(guess)
                guess = guess - np.dot(np.linalg.inv(JacobianMatrix), np.array([self.pol1([x,y]),self.pol2([x,y])])) #xn+1 = xn - J^-1*f(xn)
                xNew = guess[0]
                yNew = guess[1]
                if abs(x-xNew) <= 1.e-6 and abs(y-yNew) <= 1.e-6:   # x is still the "old" x value so if we compare it to the "new"x with subtraction and the absolute with the answer is suffieciently small (1e-6) the we conclude that the method has converged
                    if combo == True:
                        return ([xNew, yNew], i+1)
                    elif itersteps == False:
                        return [xNew, yNew]
                    else:
                        return i+1
                x = xNew
                y = yNew
        else:
            for i in range(100):
                JacobianMatrix = [np.array([self.dxPol1(guess), self.dyPol1(guess)]), np.array([self.dxPol2(guess), self.dyPol2(guess)])]       #Same as before but we do not need an approximation we just straight up put in the values of our guess in the dx, dy polinomials
                JacobianMatrix = np.array(JacobianMatrix)
                if np.linalg.det(JacobianMatrix) == 0 and abs(self.pol1(guess)) < 1.e-6 and abs(self.pol2(guess)) < 1.e-9:
                    if combo == True:
                        return (guess, i+1)
                    elif itersteps == False:
                        return guess
                    else:
                        return i+1
                elif np.linalg.det(JacobianMatrix) == 0:
                    if combo == True:
                        return (None, 0)
                    elif itersteps == False:
                        return None
                    else:
                        return 0
                elif i == 99:
                    if combo == True:
                        return (None, 0)
                    elif itersteps == False:
                        return None
                    else:
                        return 0
                JacobianMatrix = np.array(JacobianMatrix)
                guess = np.array(guess)
                guess = guess - np.dot(np.linalg.inv(JacobianMatrix), np.array([self.pol1([x,y]),self.pol2([x,y])]))
                xNew = guess[0]
                yNew = guess[1]
                if abs(x-xNew) <= 1.e-6 and abs(y-yNew) <= 1.e-6:
                    if combo == True:
                        return ([xNew, yNew], i+1)
                    elif itersteps == False:
                        return [xNew, yNew]
                    else:
                        return i + 1
                x = xNew
                y = yNew

    def simpleNewtonMethod(self, guess, itersteps = False, combo = False): # Largely the same as newtons method except the jacobian is only computed once, and outside of the main for loop. This makes it not update and makes this method extremly slow to convege
        x = guess[0]
        y = guess[1]

        if self.derivativeAprrox == True:
            JacobianMatrix = [np.array([self.approxXDer(guess), self.approxYDer(guess)]), np.array([self.approxXDer(guess), self.approxYDer(guess)])]
        else:
            JacobianMatrix = [np.array([self.dxPol1(guess), self.dyPol1(guess)]), np.array([self.dxPol2(guess), self.dyPol2(guess)])]
            
        if np.linalg.det(JacobianMatrix) == 0 and abs(self.pol1(guess)) < 1.e-6 and abs(self.pol2(guess)) < 1.e-9:
            if itersteps == False:
                return guess
            else:
                return 1
        elif np.linalg.det(JacobianMatrix) == 0:
            if itersteps == False:
                return None
            else:
                return 0
        for i in range(10000):          # Range is larger because of the horrid convergence
            if guess[0] > 1e5 or guess[1] > 1e5:    # Divergence stopper 
                if itersteps == False:
                    return None
                else:
                    return 0
            JacobianMatrix=np.array(JacobianMatrix)
            guess = np.array(guess)
            guess = guess-np.dot(np.linalg.inv(JacobianMatrix), np.array([self.pol1([x,y]),self.pol2([x,y])]))
            xNew = guess[0]
            yNew = guess[1]
            if abs(x-xNew) <= 1.e-9 and abs(y-yNew) <= 1.e-9:
                if itersteps == False:
                    return [xNew, yNew]
                else:
                    return i+1
            elif i == 9999:         
                if itersteps == False:
                    return None
                else:
                    return 0
            x = xNew
            y = yNew
        
    def NewZero(self, guess, NewtonType = 1):  #Newton type is input by the user later but all the method does os check our self.zeroes for the guess we just converged
        if NewtonType == 1:                    # IF newtontype is one use the regular newtonMethod
            Zero = self.newtonsMethod(guess)
        elif NewtonType == 2:                  # if its 2 use the simplified
            Zero = self.simpleNewtonMethod(guess)
        else:                                  # If nothing else seems to be put in, use the normal method
            Zero = self.newtonsMethod(guess)    
    

        for i in self.zeroes:                  # Loops through our list self.zeroes
            if type(Zero) == type(None):       # If the type we are checking with we knwo that the index is 0 because None is the first index
                return 0
            elif type(Zero) == list and type(i) == list:    # If the type is list (Ie we get a point [x,y]) and the type of the throwaway is also list
                if abs(Zero[0]-i[0]) <= 1.e-4 and abs(Zero[1]-i[1]) <= 1.e-4:   # Compare the zeros first point with the i[0]th position in the self.zeroes list then 
                    return self.zeroes.index(i)                                 # We have found a zero which was already in self.zeroes and then we return the index number of that zero
                elif self.zeroes.index(i)+1 == len(self.zeroes):                # else if self.zeroes.index for i (+1 to compensate for the None in the start) = the length of the entire self.zeroes (+None)
                    self.zeroes.append(Zero)                                    # We have found a new zero which we append to the list
                    return self.zeroes.index(Zero)                              # Then we return the index value for that zero
            elif len(self.zeroes) == 1:                                         # If its the first time we look through the list (and we do not have None as the zero) we know that the list just has None and nothing else
                self.zeroes.append(Zero)                                        # Since the list is empty except None we just append the zero
                return self.zeroes.index(Zero)                                  # Return the index of the Zero
            
            
            
    def Plot(self, N, a, b, c, d, UserNewtonType):             # Plot is made with 4 params. N = resolution a is the minimum x, b is the max of x, c is the min of y and d is the max of y. It also takes which newtonmethod to use
        xvalues = np.linspace(a, b, N)         # Makes a linspace of a,b witht he step size N for the x values
        yvalues = np.linspace(c, d, N)         # and y values
        [xv, yv] = np.meshgrid(xvalues, yvalues)        #Makes a meshgrid of our two lispace lists of xv and yv (stands for xvalues and yvalues)
        
        ZeroIndexMatrix = np.zeros((N,N))                                       # Creates an empty NxN matrix

        for i in range(N):                                                      
            for j in range(N):
                Point = np.array([xv[i][j], yv[i][j]]).T                        # I is the row index and J is the column index so we go trhough each point on a squraed area and takes each point
                #print(Point)       # To see if the method is working, removed if performance is paramount
                ZeroIndexMatrix[i][j] = self.NewZero(Point, UserNewtonType)     # Changes the point on the squear grid zeroindexmatrix the the number of which zero we got out of NewZero

        plt.pcolor(xv, yv, ZeroIndexMatrix)                                     # The amount of colors in pcolor is different on each point with different zero index. Goes through ZeroINdexMatix to do this
        plt.show()


        
    def NewIter(self, guess, NewtonType = 1):               # This is the same as NewZero but instead of the point we are lookin for the number of iterations it took to converge our point
            
            if NewtonType == 1:
                Iter = self.newtonsMethod(guess,True)       # Sends the guess but the itersteps argument is true making the method return a number instead of a zero
            elif NewtonType == 2:
                Iter = self.simpleNewtonMethod(guess, True)
            else:
                Iter = self.newtonsMethod(guess,True)

            if Iter in self.zeroIndex:                      # the checks doesnt need to be as complicated as we are just comparing a number to another each time (None is just zero after all)
                return self.zeroIndex.index(Iter)
            else:
                self.zeroIndex.append(Iter)


    def NewCombo(self, guess, NewtonType = 1):  #Newton type is input by the user later but all the method does os check our self.zeroes for the guess we just converged
        if NewtonType == 1:                    # IF newtontype is one use the regular newtonMethod
            Zero, Iter = self.newtonsMethod(guess, combo=True)
        elif NewtonType == 2:                  # if its 2 use the simplified
            Zero, Iter = self.simpleNewtonMethod(guess, combo=True)
        else:                                  # If nothing else seems to be put in, use the normal method
            Zero, Iter = self.newtonsMethod(guess, combo=True)    


        iterIndex = 0
        
        if Iter in self.zeroIndex:                      # the checks doesnt need to be as complicated as we are just comparing a number to another each time (None is just zero after all)
            iterIndex = self.zeroIndex.index(Iter)
        else:
            self.zeroIndex.append(Iter)
            iterIndex = self.zeroIndex.index(Iter)
                        
        zeroIndex = 0
        
        for i in self.zeroes:                  # Loops through our list self.zeroes
            if type(Zero) == type(None):       # If the type we are checking with we knwo that the index is 0 because None is the first index
                zeroIndex = 0
            elif type(Zero) == list and type(i) == list:    # If the type is list (Ie we get a point [x,y]) and the type of the throwaway is also list
                if abs(Zero[0]-i[0]) <= 1.e-4 and abs(Zero[1]-i[1]) <= 1.e-4:   # Compare the zeros first point with the i[0]th position in the self.zeroes list then 
                    zeroIndex = self.zeroes.index(i)                                 # We have found a zero which was already in self.zeroes and then we return the index number of that zero
                elif self.zeroes.index(i)+1 == len(self.zeroes):                # else if self.zeroes.index for i (+1 to compensate for the None in the start) = the length of the entire self.zeroes (+None)
                    self.zeroes.append(Zero)                                    # We have found a new zero which we append to the list
                    zeroIndex = self.zeroes.index(Zero)                              # Then we return the index value for that zero
            elif len(self.zeroes) == 1:                                         # If its the first time we look through the list (and we do not have None as the zero) we know that the list just has None and nothing else
                self.zeroes.append(Zero)                                        # Since the list is empty except None we just append the zero
                zeroIndex = self.zeroes.index(Zero)                                  # Return the index of the Zero
            
        return (zeroIndex, iterIndex)

                    
            
            

    def plotItr(self, N, a, b, c, d, UserNewtonType): # Only difference between this plot and the normal Plot is that we can get a different color for each zero and how long(many iterations) it took to converge
        xvalues = np.linspace(a, b, N)
        yvalues = np.linspace(c, d, N)
        [xv, yv] = np.meshgrid(xvalues, yvalues)

        IterIndexMatrix = np.zeros((N,N))
        
        for i in range(N):
            for j in range(N):
                Point = np.array([xv[i][j], yv[i][j]]).T
                #print(self.NewIter(Point, UserNewtonType))
                IterIndexMatrix[i][j] = self.NewIter(Point, UserNewtonType)

        #print(IterIndexMatrix)
        plt.pcolor(xv, yv, IterIndexMatrix)
        plt.show()



    def PlotCombo(self, N, a, b, c, d, UserNewtonType=1): 
        xvalues = np.linspace(a, b, N)
        yvalues = np.linspace(c, d, N)
        [xv, yv] = np.meshgrid(xvalues, yvalues)

        IterIndexMatrix = np.zeros((N,N))
        ZeroIndexMatrix = np.zeros((N,N)) 
        
        for i in range(N):
            for j in range(N):
                Point = np.array([xv[i][j], yv[i][j]]).T
                IterIndexMatrix[i][j], ZeroIndexMatrix[i][j] = self.NewCombo(Point, UserNewtonType)

        #print(IterIndexMatrix)
        my_cmap = cm.gray
        plt.pcolormesh(xv, yv, ZeroIndexMatrix)
        plt.pcolormesh(xv, yv, IterIndexMatrix, cmap=my_cmap, alpha=0.3)
        plt.show()
                
        
# The three combinations of our functions:
p = fractal2D(F, G, dxF, dyF, dxG, dyG) 
p2 = fractal2D(F, G) 
q = fractal2D(H, I, dxH, dyH, dxI, dyI)
q2 = fractal2D(H, I)
r = fractal2D(J, K, dxJ, dyJ, dxK, dyK)
r2 = fractal2D(J, K)


#Change p to q or r for different polinomials
r.Plot(100, -7.5, 7.5, -7.5, 7.5, 1)
#r.PlotCombo(100, -7.5, 7.5, -7.5, 7.5, 1)

