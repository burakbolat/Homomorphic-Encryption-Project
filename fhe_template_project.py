from re import I
import re
from eva import EvaProgram, Input, Output, evaluate, py_to_eva
from eva.ckks import CKKSCompiler
from eva.seal import generate_keys
from eva.metric import valuation_mse
import timeit
import networkx as nx
from random import Random

import matplotlib.pyplot as plt
import numpy as np

N = -1  # Global number of nodes for program analitic program to know number of nodes

def DFS(v, visited, graph):
    visited[v] = True
    for neighbour in graph[v]:
        if not visited[neighbour]:
            DFS(neighbour, visited, graph)

# Using networkx, generate a random graph
# You can change the way you generate the graph
def generateGraph(n, p):
    G = nx.erdos_renyi_graph(n, p)
    # G = nx.Graph()
    # G.add_nodes_from([0, 1, 2, 3, 4, 5])
    # G.add_edges_from([(0,1), (0,2), (1,3), (2,3), (4,5)])
    
    nx.draw(G, with_labels=True)
    plt.savefig('graph.png')
    return G

# If there is an edge between two vertices its weight is 1 otherwise it is zero
# You can change the weight assignment as required
# Two dimensional adjacency matrix is represented as a vector
# Assume there are n vertices
# (i,j)th element of the adjacency matrix corresponds to (i*n + j)th element in the vector representations
def serializeGraphZeroOne(GG,vec_size):
    n = GG.size() # Returns the number of edges or total of all edge weights.
    n = len(GG.nodes)
    graphdict = {}
    g = []
    for row in range(n):
        for column in range(n):
            if GG.has_edge(row, column):
                weight = 1.0
            else:
                weight = 0.0
            g.append( weight  )  
            key = str(row)+'-'+str(column)
            graphdict[key] = [weight] # EVA requires str:listoffloat
    # EVA vector size has to be large, if the vector representation of the graph is smaller, fill the eva vector with zeros
    for i in range(vec_size - n*n): 
        g.append(0.0)
    return g, graphdict

# To display the generated graph
def printGraph(graph,n):
    for row in range(n):
        for column in range(n):
            print("{:.5f}".format(abs(graph[row*n+column])), end = '\t')  # Added abs since - sign breaks visible continuity
        print() 

# Eva requires special input, this function prepares the eva input
# Eva will then encrypt them
def prepareInput(n, m):
    input = {}
    GG = generateGraph(n, 0.07)

    visited = [False] * len(GG.nodes)
    head = 0
    DFS(head, visited, GG)

    print("Ground Truth:", list(nx.connected_components(GG)))
    # for ind in range(len(visited)):
    #     if not visited[ind]: print(ind)

    graph, graphdict = serializeGraphZeroOne(GG,m)
    input['Graph'] = graph
    return input

def get_mask(graph, j, n):
    y = []
    for k in range(4096*4):
        if k == j: y.append(1.0)
        else: y.append(0.0)
    y = py_to_eva(y)
    y = y*graph

    i = 1
    while i < n/2:
        y += y>>i
        i <<= 1
        
    if i < n:
        y += y>>(n-i)
        
    y = y>>(n*j-j)
    return y

def get_filtered_neigs(graph, N):
    masks = []
    for j in range(N):
        mask = get_mask(graph, j, N)
        masks.append(mask)

    sum_of_masks = graph.program.vec_size * [0.0]
    for mask in masks:
        sum_of_masks += mask
    filtered_neighbours = sum_of_masks*graph
    return filtered_neighbours

## ! Not used in current implementation
## Planned for matrix multiplication implementation
def row_times_col(mat1, mat2, row, col):
    global N

    matmul = (mat1<<(row*N)) * (mat2<<(col*N))
    temp = []
    while len(temp) < 4096*4:
        temp.append(0.0)
    for i in range(N):
        temp += (matmul << (N-i))
    temp = temp >> ((row*N) + col)
    return temp

def graphanalticprogram(graph):
    global N
    
    # for m in range(4):
    filtered_neighbours = get_filtered_neigs(graph, N)
    for k in range(1,N):
        graph += (filtered_neighbours << (N*k))

    return graph

    # reval = graph
    # # while len(reval) < 4096*4:
    # #     reval.append(0.0)
    # identity = []
    # for i in range(4):
    #     for j in range(4):
    #         if i==j: identity.append(1.0)
    #         else: identity.append(0.0)
    # while len(identity) < 4096*4:
    #     identity.append(0.0)

    # n=6
    # # for deg in range(n): 
    # for row in range(1):
    #     for col in range(n):
    #         temp = row_times_col(reval, graph, row, col)
    #         mask = []
    #         for mask_ind in range(4096*4):
    #             if mask_ind == row*n + col : mask.append(1.0)
    #             else: mask.append(0.0)
    #         reval += temp*mask

    # for i in range(4):
    #     reval += identity
    #     reval = reval *
     
   
    # return reval

# Do not change this 
# the parameter n can be passed in the call from simulate function
class EvaProgramDriver(EvaProgram):
    def __init__(self, name, vec_size=4096, n=4):
        self.n = n
        super().__init__(name, vec_size)

    def __enter__(self):
        super().__enter__()

    def __exit__(self, exc_type, exc_value, traceback):
        super().__exit__(exc_type, exc_value, traceback)

# Repeat the experiments and show averages with confidence intervals
# You can modify the input parameters
# n is the number of nodes in your graph
# If you require additional parameters, add them
def simulate(n):
    m = 4096*4
    global N 
    N = n
    print("Will start simulation for ", n)
    config = {}
    config['warn_vec_size'] = 'false'
    config['lazy_relinearize'] = 'true'
    config['rescaler'] = 'always'
    config['balance_reductions'] = 'true'
    inputs = prepareInput(n, m)
    graphanaltic = EvaProgramDriver("graphanaltic", vec_size=m,n=n)
    old_head = inputs["Graph"][:N]
    for _ in range(N-1):

        with graphanaltic:
            graph = Input('Graph')
            reval = graphanalticprogram(graph)
            Output('ReturnedValue', reval)
        
        prog = graphanaltic
        prog.set_output_ranges(30)
        prog.set_input_scales(60)

        start = timeit.default_timer()
        compiler = CKKSCompiler(config=config)
        compiled_multfunc, params, signature = compiler.compile(prog)
        compiletime = (timeit.default_timer() - start) * 1000.0 #ms

        start = timeit.default_timer()
        public_ctx, secret_ctx = generate_keys(params)
        keygenerationtime = (timeit.default_timer() - start) * 1000.0 #ms
        
        start = timeit.default_timer()
        encInputs = public_ctx.encrypt(inputs, signature)
        encryptiontime = (timeit.default_timer() - start) * 1000.0 #ms

        start = timeit.default_timer()
        encOutputs = public_ctx.execute(compiled_multfunc, encInputs)
        executiontime = (timeit.default_timer() - start) * 1000.0 #ms

        start = timeit.default_timer()
        outputs = secret_ctx.decrypt(encOutputs, signature)
        decryptiontime = (timeit.default_timer() - start) * 1000.0 #ms

        start = timeit.default_timer()
        reference = evaluate(compiled_multfunc, inputs)
        referenceexecutiontime = (timeit.default_timer() - start) * 1000.0 #ms
        ret_val = np.array(outputs["ReturnedValue"])
        ret_val = [1 if ret_val[i]>0.5 else 0 for i in range(ret_val.shape[0])]
        if(np.allclose(old_head, ret_val[:N])): break
        else: old_head = ret_val[:N] 
        # print("modified", ret_val[:N])
        inputs["Graph"] = ret_val
    print("Executed:", ret_val[:N])
    print("Connected", np.count_nonzero(ret_val[:N]), "Non-connected:", N-np.count_nonzero(ret_val[:N]))
    
    # Change this if you want to output something or comment out the two lines below
    # for key in outputs:
        # print(key, float(outputs[key][0]), float(reference[key][0]))
        # printGraph(outputs[key], n)
        # output = outputs[key][:N*N]
        # output = np.array(output[:N])
        # print(output)
        # printGraph(outputs[key], n)
        # print("----------------")
        # printGraph(inputs['Graph'], n)
        # print("----------------")
        # printGraph(reference[key], n)

    mse = valuation_mse(outputs, reference)  # since CKKS does approximate computations, this is an important measure that depicts the amount of error

    return compiletime, keygenerationtime, encryptiontime, executiontime, decryptiontime, referenceexecutiontime, mse


if __name__ == "__main__":
    simcnt = 1  # The number of simulation runs, set it to 3 during development otherwise you will wait for a long time
    # For benchmarking you must set it to a large number, e.g., 100
    # Note that file is opened in append mode, previous results will be kept in the file
    # resultfile = open("results.csv", "a")  # Measurement results are collated in this file for you to plot later on
    # resultfile.write("NodeCount,SimCnt,CompileTime,KeyGenerationTime,EncryptionTime,ExecutionTime,DecryptionTime,ReferenceExecutionTime,Mse\n")
    # resultfile.close()
    
    print("Simulation campaing started:")
    for nc in range(36,64,4): # Node counts for experimenting various graph sizes
        resultfile = open("results.csv", "w") 
        for i in range(simcnt):
            # Call the simulator
            compiletime, keygenerationtime, encryptiontime, executiontime, decryptiontime, referenceexecutiontime, mse = simulate(nc)
            res = str(nc) + "," + str(i) + "," + str(compiletime) + "," + str(keygenerationtime) + "," +  str(encryptiontime) + "," +  str(executiontime) + "," +  str(decryptiontime) + "," +  str(referenceexecutiontime) + "," +  str(mse) + "\n"
            
            ## For milestone 4

            # str1 = ""
            # for i in range(n):
            #     for j in range(n):
            #         # val = 1 if output[i*n+j] > 0.8 else 0
            #         val = output[i*n+j]
            #         str1 += f'{val:.3f}' + " "
            #     str1 += "\n" 
            # resultfile.write(str1)
            # resultfile.close()
