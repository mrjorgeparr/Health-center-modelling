#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 20:49:06 2020

@author: isegura
"""
import sys


class HealthCenter():
    def __init__(self,name=None):
        self.name=name
        
        
    def __eq__(self,other):
        return  other!=None and self.name==other.name
    
    def __str__(self):
        return self.name


class AdjacentVertex:
  def __init__(self,vertex,weight):
    self.vertex=vertex
    self.weight=weight
  
  def __str__(self):
    return '('+str(self.vertex)+','+str(self.weight)+')'
 
class Map():
    def __init__(self):
        self.centers={}
        self.vertices={}
    
    def addHealthCenter(self,center):
        i=len(self.centers)
        self.centers[i]=center
        self.vertices[i]=[]
        
    def _getIndice(self,center):
        for index in self.centers.keys():
            if self.centers[index]==center:
                return index
        return -1
        
    def __str__(self):
        result=''
        for i in self.vertices.keys():
            result+=str(self.centers[i])+':\n'
            for adj in self.vertices[i]:
                result+='\t'+str(self.centers[adj.vertex])+', distance:'+str(adj.weight)+'\n'
        return result

    def areConnectedv2(self, center1, center2, distance):
        if type(center1) is not HealthCenter or type(center2) is not HealthCenter:
            print("Error in input variable/s!")
            return
        index1 = self._getIndice(center1)
        index2 = self._getIndice(center2)
        if index1 == -1:
            print(center1, ' not found!')
            return
        if index2 == -1:
            print(center2, ' not found!!')
            return

        for adj in self.vertices[index1]:
            if adj.vertex == index2 and adj.weight == distance:
                return adj.weight
        # print(pto1,pto2," no están conectados")
        return 0

    def addConnection(self,center1,center2,distance):
        #print('new conexion:',pto1,pto2)
        index1=self._getIndice(center1)
        index2=self._getIndice(center2)
        if index1==-1:
            print(center1,' not found!')
            return
        if index2==-1:
            print(center2,' not found!!')
            return
            #we must check if there's already a connection of that distance
        self.vertices[index1].append(AdjacentVertex(index2,distance))
        #print('adding:',index2,index1,distancia)
        self.vertices[index2].append(AdjacentVertex(index1,distance))


        
    def areConnected(self,center1,center2):
        index1=self._getIndice(center1)
        index2=self._getIndice(center2)
        if index1==-1:
            print(center1,' not found!')
            return
        if index2==-1:
            print(center2,' not found!!')
            return 
        
        for adj in self.vertices[index1]:
            if adj.vertex==index2:
                return adj.weight
        #print(pto1,pto2," no están conectados")
        return 0
            
    def removeConnection(self,center1,center2):
        index1=self._getIndice(center1)
        index2=self._getIndice(center2)
        if index1==-1:
            print(center1,' not found!')
            return
        if index2==-1:
            print(center2,' not found!!')
            return 
        
        for adj in self.vertices[index1]:
            if adj.vertex==index2:
                self.vertices[index1].remove(adj)
                break
                
        for adj in self.vertices[index2]:
            if adj.vertex==index1:
                self.vertices[index2].remove(adj)
                break

    def removeConnection_weight(self, center1, center2, weight):
        index1 = self._getIndice(center1)
        index2 = self._getIndice(center2)
        if index1 == -1:
            print(center1, ' not found!')
            return
        if index2 == -1:
            print(center2, ' not found!!')
            return

        for adj in self.vertices[index1]:
            if adj.vertex == index2 and adj.weight == weight:
                self.vertices[index1].remove(adj)
                break

        for adj in self.vertices[index2]:
            if adj.vertex == index1 and adj.weight == weight:
                self.vertices[index2].remove(adj)
                break
        

    def createPath(self): 
        """This function prints the vertices by dfs algorithm"""
        #print('dfs traversal:')
        # Mark all the vertices as not visited 
        visited = [False] * len(self.vertices)

        paths=[]
        for v in  self.vertices:
            if visited[v]==False:
                self._dfs(v, visited,paths)
        
        print()
        return paths
        
    def _dfs(self, v, visited,paths): 
        # Mark the current node as visited and print it 
        visited[v] = True
        #print(self.centers[v], end = ' ') 
        paths.append(self.centers[v])
        # Recur for all the vertices  adjacent to this vertex 
        for adj in self.vertices[v]: 
          i=adj.vertex
          if visited[i] == False: 
            self._dfs(i, visited,paths) 
            
            
    
    def printSolution(self,distances,previous,v): 
        """imprime los caminos mínimos desde v"""
        for i in range(len(self.vertices)):
          if distances[i]==sys.maxsize:
            print("There is not path from ",v,' to ',i)
          else: 
            minimum_path=[]
            prev=previous[i]
            while prev!=-1:
              minimum_path.insert(0,self.centers[prev])
              prev=previous[prev]
            
            minimum_path.append(self.centers[i])  
    
            print('Ruta mínima de:',self.centers[v],'->',self.centers[i],", distance", distances[i], ', ruta: ',  end= ' ')
            for x in minimum_path:
                print(x,end= ' ')
            print()
    
    def minDistance(self, distances, visited): 
        """This functions returns the vertex (index) with the mininum distance. To do this, 
        we see in the list distances. We 
        only consider the set of vertices that have not been visited"""
        # Initilaize minimum distance for next node 
        min = sys.maxsize 
    
        #returns the vertex with minimum distance from the non-visited vertices
        for i in range(len(self.vertices)): 
          if distances[i] <= min and visited[i] == False: 
            min = distances[i] 
            min_index = i 
      
        return min_index 
    
    def dijkstra(self, v=0): 
        """"This function takes the index of a delivery point pto and calculates its mininum path 
        to the rest of vertices by using the Dijkstra algoritm. Devuelve una lista con las distancias
        y una lista con los vértices anteriores a uno dado en el camino mínimo"""  
        
        
        #we use a Python list of boolean to save those nodes that have already been visited  
        visited = [False] * len(self.vertices) 
    
        #this list will save the previous vertex 
        previous=[-1]*len(self.vertices) 
    
        #This array will save the accumulate distance from v to each node
        distances = [sys.maxsize]*len(self.vertices) 
        #The distance from v to itself is 0
        distances[v] = 0
    
        for i in range(len(self.vertices)): 
          # Pick the vertex with the minimum distance vertex.
          # u is always equal to v in first iteration 
          u = self.minDistance(distances, visited) 
          # Put the minimum distance vertex in the shotest path tree
          visited[u] = True
          
          # Update distance value of the u's adjacent vertices only if the current  
          # distance is greater than new distance and the vertex in not in the shotest path tree 
          for adj in self.vertices[u]:
            i=adj.vertex
            w=adj.weight
            if visited[i]==False and distances[i]>distances[u]+w:
              distances[i]=distances[u]+w   
              previous[i]=u       
              
        #finally, we print the minimum path from v to the other vertices
        #self.printSolution(distances,previous,v)
        return previous,distances
 
    def minimumPath(self, start, end):
        """calcula la ruta mínima entre dos puntos de entrega"""
        indexStart=self._getIndice(start)
        if indexStart==-1:
            print(str(start) + " does not exist")
            return None
        indexEnd=self._getIndice(end)
        if indexEnd==-1:
            print(str(end)  + " does not exist")
            return None
        
        previous,distances=self.dijkstra(indexStart)

        #construimos el camino mínimo
        minimum_path=[]
        prev=previous[indexEnd]
        while prev!=-1:
            minimum_path.insert(0,self.centers[prev])
            prev=previous[prev]
            
        minimum_path.append(self.centers[indexEnd])
        return minimum_path, distances[indexEnd]

    def error_checking(self, start, end):
        if (type(start) is not HealthCenter) or (type(end) is not HealthCenter):
            print("\nError: input variables were of invalid type")
            return -1,-1
        index = self._getIndice(start)
        index2 = self._getIndice(end)
        if index == -1:
            print("\nStarting node was not on the map")
            return -1,-1
        if index2 == -1:
            print("\nGoal node was not on the map")
            return -1, -1
        return index, index2

    def minimumPathBF(self , start, end):
        """Regarding complexity: O(|V|*|E|) <---- worst case performance, where |V| := number of vertices, and 
        |E| := number of edges, expressed in terms of  the numbers of edges and vertices.
        The function is clearly O(n^3) as in the main part of Bellman-ford, by definition we must check a condition for 
        each vertex, and for each connection that the vertex has, which implies a nested loop, and thus the complexity 
        is n^3, since we check until convergence, this algorithm is a bit slower than Dijkstra but it can handle 
        negative weights, which is not our case since distances are always positive. 
        We also find nested loops in the part relating to checking to see whether the algorithm converges or not.
        
        The first thing we do is check if the input variables are of the correct type, with a simple if statement, then,
        assuming correct input we search for the health centers, to see if they are on the invoking Map or not. Those 
        are the best cases, as they are constant time operations.
        
        The worst cases are when the function reaches the checking part, and the algorithm converges, in which the 
        function does two nested loops and a linear time loop, only when it converges, to restore the path.
        """
        minimum_path = []
        minimumDistance = 0
        index, index2 = self.error_checking(start, end)
        if index == -1 and index2 == -1:
            return minimum_path, minimumDistance
        #for each vertex we initialize d and p
        d = [float('inf')]*len(self.centers)
        p = [None]*len(self.centers)
        p[index] = 0
        d[index] = 0
        for j in range(len(self.centers)):
            for i in range(len(self.centers)):          #for each vertex
                if i == index:
                    continue
                ls = self.vertices[i]            # we consider the list of adjacency
                v = i
                for con in ls:                   #and iterate through it
                    u = con.vertex
                    w = con.weight
                    if w < 0:       # we do not allow negative weights
                        continue
                    if d[v] > (d[u] + w):
                        d[v] = d[u] + w
                        p[v] = u
            # we must check all connections
        for i in range(len(self.centers)):          #for each vertex
            if i == index:
                continue
            ls = self.vertices[i]            # we consider the list of adjacency
            v = i
            for con in ls:                   #and iterate through it
                u = con.vertex
                w = con.weight
                if d[v] > (d[u] + w):   #u is the index for the health center considered at every iteration
                    # the algorithm does not converge
                    print("does not converge")
                    return False,-1
        n = index2
        minimumDistance = d[n]       #now we need to restore the path from p and d
        while n != index:
            minimum_path.append(str(self.centers[p[n]]))
            n = p[n]
        minimum_path.reverse()      #we must reverse it as we were storing them backwards on the loop
        minimum_path.append(str(self.centers[index2]))
        return minimum_path, minimumDistance

    def minimumPathFW(self, start , end):
        """Floyd-Warshall (G)
            Initialize
            Let D be the adjacency matrix (distances between vertices)
            Let P be a matrix to store the minimum paths
            If i=j or Dij= infinite then: 
                Pi,j = null
            otherwise: 
                Pi,j=i 
                
            For k = 0 To len(V)-1
                For i = 0 To len(V)-1
                    if i==k: continue
                    For j = 0 To len(V)-1
                        Di,j = min(Di,j , Di,k + Dk,j )
                        If min = Di,k + Dk,j then
                        Pi,j = Pk,j                      
        """
        """
        Regarding complexity:
        The first thing we do is check for invalid inputs, maybe invalid type or non existing health centers, then 
        assuming no bad inputs, we declare the matrices d, and p, with every element initialized as infinity and 
        None respectively.
        
        We fill the matrices with the appropriate values for the distances and the paths. In the matrix D, we must place 
        all the values from all connections , and in the diagonals set to 0, as there are no loops on map, for the 
        matrix P, the values in the diagonals are set to None, as we are not interested in the path from a loop to 
        itself, because as we mentioned before, there are no loops. The complexity of this part is O(n^2).
        
        Then we enter the main loop, in which we find three loops, one inside the other, it is important to skip  all 
        the iterations in which we find that j == k, as that results in unwanted changes in the p-matrix, because in 
        this case d[i][j], d[i][k] + d[j][k] since d[j][k] = d[j][j] = 0 and d[i][j] = d[i][k], so they are equal but 
        there is no change in the d-matrix, and thus there should be no change in the p-matrix either, if we did not 
        take that into account, the entire p-matrix would result in None, also, in the inner loop, we remove all the 
        iterations in which i == k, to avoid redundant comparisons.
        
        In this way we reduce the number of iterations (a total of |V|^3, where |V| is the number of vertices), 
        considerably as instead, we find a total of |V|*(|V|-1)*(|V|-2), which means a reduction of 2*|V|*(2*|V| - 1) 
        iterations. 
        
        The last part is the reconstruction of the path, from the values on the matrix P. This part is linear, as it 
        iterates through all the values in the row corresponding to the start input parameter (s), beginning at the position 
        corresponding to the end input parameter (e), until we reach column corresponding to s.
        
                
        In terms of complexity we find that the fastest case scenario, is when the input is not valid, in which the 
        complexity is constant, and the worst case scenario is when the algorithm reaches its main part where we can 
        find three loops, which gives us a complexity of O(n^3), although there's some reduction as some iterations are 
        skipped. 
        """
        minimum_path=[]
        minimumDistance=0
        index, index2 = self.error_checking(start, end)         # we check for invalid input, that is, possible syntactic errors, or non existant health centers
        if index == -1 and index2 == -1:
            return minimum_path, minimumDistance
        d = [[float('inf') for col in range(len(self.centers))] for row in range(len(self.centers))]            #we declare the matrices as required
        p = [[None for col in range(len(self.centers))] for row in range(len(self.centers))]

        for v in range(len(self.centers)):          #we fill the matrix D, with the values for the distances
            d[v][v] = 0                     #all the diagonal elements have value 0 in D, as D(i,j) = 0
            ls = self.vertices[v]
            for con in ls:
                u = con.vertex
                d[v][u] = con.weight
                p[v][u] = v

        for k in range(len(self.centers)):
            for j in range(len(self.centers)):
                if j == k:                    # if j == k, d(j,k) is 0, and d(i,j) == d(i,k), which leads to error
                    continue
                for i in range(len(self.centers)):
                    if i == k or i == j:  # to remove redundant changes in the p-matrix
                        continue
                    minim = d[i][j]
                    d[i][j] = min(d[i][j], d[i][k] + d[j][k])
                    if d[i][j] == d[i][k] + d[j][k] and d[i][j] != minim:
                        p[i][j] = p[k][j]

        #now we reconstruct the path
        n = index2
        while n != index:
            minimum_path.insert(0, str(self.centers[p[index][n]]))
            n = p[index][n]
        minimum_path.append(str(self.centers[index2]))
        minimumDistance = d[index][index2]
        return minimum_path, minimumDistance

def test():
    #https://www.bogotobogo.com/python/images/Dijkstra/graph_diagram.png
   
    m=Map()
    for c in ['A','B','C','D','E','F']:
        m.addHealthCenter(HealthCenter(c))
    
    print(m)

    m.addConnection(m.centers[0],m.centers[1],7)#A,B,7
    m.addConnection(m.centers[0],m.centers[2],9)#A,C,9
    m.addConnection(m.centers[0],m.centers[5],14)#A,F,14
    m.addConnection(m.centers[0], m.centers[0], -10)
    m.addConnection(m.centers[1],m.centers[2],10)#B,C,10
    m.addConnection(m.centers[1],m.centers[3],15)#B,D,15
    
    m.addConnection(m.centers[2],m.centers[3],11)#C,D,11
    m.addConnection(m.centers[2],m.centers[5],2)#C,F,2
    
    m.addConnection(m.centers[3],m.centers[4],6)#D,E,6
    
    m.addConnection(m.centers[4],m.centers[5],9)#E,F,9
    print(m)

    
    c1=m.centers[0]
    c2=m.centers[3]
    print(c1,c2,' are connected?:',m.areConnected(c1,c2))
    
    c2=m.centers[1]
    print(c1,c2,' are connected?:',m.areConnected(c1,c2))

    m.removeConnection(c1,c2)
    print(m)
    
    print('createPath:',end=' ')
    ruta=m.createPath()
    #print('Ruta:',ruta)
    for r in ruta:
        print(r, end=' ')
    print()

    minimum_path,d=m.minimumPath(c1,m.centers[5])
    for p in minimum_path:
        print(p,end=' ')
    print('total distance:',d)

    #añade más pruebas para probar los dos nuevos métodos minimumPathBF y minimumPathFW
    minimum_path,d=m.minimumPathBF(m.centers[0],m.centers[0])

    for p in minimum_path:
        print(p,end=' ')
    print('total distance:',d)
    minimum_path,d=m.minimumPathFW(c1, c2)
    print("Distancia: ", d)
    print("Camino: ",minimum_path)
    """
    for p in minimum_path:
        print(p,end=' ')
    print('total distance:',d)
    """


#Descomenar para usarlo en Spyder
test()
