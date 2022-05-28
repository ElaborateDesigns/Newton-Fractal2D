import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import matplotlib.colors as col
import sys



###################
### CLASS START ###
###################

    
class fractal2D:
    def __init__(self, pol1, pol2, nrIterations = 10, err = 1e-100, 
                 min_x = 0,
                 min_y = 0,
                 max_x = 3,
                 max_y = 3,
                 size_x = 4, 
                 size_y = 3):
        self.pol1 = pol1
        self.pol2 = pol2
        self.nrIterations = nrIterations
        self.err = err

        if pol1 == pol2:
            raise SyntaxError("You have entered the same function")
        
        #Grid Calculations stuff             
        x = np.linspace (min_x, max_x, size_x)
        y = np.linspace (min_y, max_y, size_y)
        self.grid = np.meshgrid(x,y)
        print(self.grid)
        self.A = np.empty( (size_x, size_y) )
        self.B = np.zeros(self.A.shape) #default


        
#######################
### CALCULATE ROOTS ###
###     SECTION     ###
#######################        


    def newtonsMethod(self, guess):
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
        
        #print(jacobianMatrix)
        for i in range(2): #loop to lambdify the functions in the jacobianMatrix(from sympy to function)
            jacobianMatrix[i,0] = sym.lambdify([a,b], jacobianMatrix[i,0], 'numpy')
            jacobianMatrix[i,1] = sym.lambdify([a,b], jacobianMatrix[i,1], 'numpy')
        #print(jacobianMatrix[1,1](1,1))
        
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
            if np.linalg.det(jMnumbers) == 0:
                jMnumbers[0,0] = jMnumbers[0,0] + self.err
                jMnumbers[1,1] = jMnumbers[0,1] + self.err
            jMnumbers = np.linalg.inv(jMnumbers)
            
            dotProduct = np.dot(jMnumbers,funMatrix)
            guess = guess - dotProduct #xn-1 = xn - J^-1*f(xn)
            gList.append(guess)
            if len(gList) > 2:
                if abs(gList[i][0] - gList[i-1][0]) < self.err and abs(gList[i][1] - gList[i-1][1]) < self.err:
                    #print(guess)
                    return guess

            #print(i)
        

        #print(guess)
        return guess
    
    
    
    
    def getAmountZeros(self,x0):
            
        Zeros=[]
        Tol = 1e-16
        print(len(x0))
        if len(x0) > 1:
            for i in range(len(x0)):
                x1 = [x0[i][0], x0[i][1]]
                print(i)
                print(x1)
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
        
        
    
    def plot(self):
        xmin, ymin = -10, -10
        xmax, ymax = 10, 10
        nx, ny = (4, 3)
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        xv, yv = np.meshgrid(x, y)
        self.getAmountZeros(yv)
        


    
    def plot2(self): #This is bad and inefficient but I just wanted a working prototype. 
        xmin, ymin = -10, -10
        xmax, ymax = 10, 10
        nx, ny = (50, 50)
        x = np.linspace(xmin, xmax, nx)
        y = np.linspace(ymin, ymax, ny)
        xv, yv = np.meshgrid(x, y)
        #self.getAmountZeros(yv)
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
                result = self.newtonsMethod( (x0,y0) )
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
      palette = [(2,0.5,0.5), (0.5,2.0,0.5), (0.5,0.5,2), (3,3,0), (0, 2, 2), (2, 0, 2) ]
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
p = fractal2D(a**8 - 28*a**6*b**2 + 70*a**4 + 15*a**4 - 28*a**2*b**6 - 90*a**2 + b**8 + 15*b**4 - 16, 8*a**7*b -56*a**5*b**3 + 56*a**3*b**5 + 60*a**3*b - 8*a*b**7 - 60*a*b**3, 1000)

###################################################################################Polies#####################################################################################

#a**8 - 28*a**6*b**2 + 70*a**4 + 15*a**4 - 28*a**2*b**6 - 90*a**2 + b**8 + 15*b**4 - 16
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

fractal2D.plot2(p)
#fractal2D.getAmountZeros(p,yv)