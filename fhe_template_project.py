from eva import EvaProgram, Input, Output, evaluate, py_to_eva
from eva.ckks import CKKSCompiler
from eva.seal import generate_keys
from eva.metric import valuation_mse
import timeit
import networkx as nx

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
def serializeGraphZeroOne(GG,vec_size, for_matmul=False):
    n = GG.size() # Returns the number of edges or total of all edge weights.
    n = len(GG.nodes)
    graphdict = {}
    g = []
    
    if for_matmul:
        for row in range(n):
            for column in range(n):
                if GG.has_edge(row, column):
                    weight = 1.0
                else:
                    weight = 0.0
                g.append(weight)  
                key = str(row)+'-'+str(column)
                graphdict[key] = [weight] # EVA requires str:listoffloat
        
    for row in range(n):
        for column in range(n):
            if GG.has_edge(row, column):
                weight = 1.0
            else:
                weight = 0.0
            g.append(weight)  
            key = str(row)+'-'+str(column)
            graphdict[key] = [weight] # EVA requires str:listoffloat
    # EVA vector size has to be large, if the vector representation of the graph is smaller, fill the eva vector with zeros
    for i in range(vec_size - len(g)): 
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
def prepareInput(n, m, p=0.3, for_matmul=False):
    input = {}
    GG = generateGraph(n, p)

    # visited = [False] * len(GG.nodes)
    # head = 0
    # DFS(head, visited, GG)

    print("Ground Truth:", list(nx.connected_components(GG)))
    print("num-of-con:", len(list(nx.connected_components(GG))[0]))
    # for ind in range(len(visited)):
    #     if not visited[ind]: print(ind)

    graph, graphdict = serializeGraphZeroOne(GG,m, for_matmul)
    input['Graph'] = graph
    return input

### Graph Analitc Program functions
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

def graphanalticprogram(graph):
    global N
    
    filtered_neighbours = get_filtered_neigs(graph, N)
    for k in range(1,N):
        graph += (filtered_neighbours << (N*k))

    return graph

# Belov is an implementation of
# Efficient matrix multiplication
# https://eprint.iacr.org/2018/1041.pdf
# https://github.com/microsoft/SEAL/issues/137
def get_u_sigma(k, d): 
    eps = 1e-6  # Read https://github.com/microsoft/SEAL/issues/137 
    n = d*d

    u_sigma = []
    if k>=0:    
        for l in range(n):
            u_sigma.append(1.0) if ((l-(d*k) < (d-k)) and (l-(d*k) >=0)) else u_sigma.append(eps)
    else:
        for l in range(n):
            u_sigma.append(1.0) if (-k <= (l-(d+k)*d) and (l-(d+k)*d) < d) else u_sigma.append(eps)

    while len(u_sigma) < 4096*4:
        u_sigma.append(0)

    return u_sigma

def get_u_tau(k, d):
    n = d*d
    u_tau = [0]*4096*4
    for i in range(d):
        u_tau[k+d*i] = 1.0

    return u_tau


def get_cta(ct, d):
    ctA = ct*get_u_sigma(0, d)

    for k in range(-d+1, d):
        if k == 0: continue  # Computed just before the loop      
        u_sigma = get_u_sigma(k, d)
        rot_val = k if k > 0 else d*d+k
        ctA = ctA + (ct<<rot_val)*u_sigma

    return ctA
        
def get_ctb(ct, d):
    ctB = ct*get_u_tau(0, d)
    for k in range(1, d):
        u_tau = get_u_tau(k, d)
        ctB += (ct<<(k*d))*u_tau

    return ctB

def get_ak(ct, k, d):
    vk = [1e-6]*(4096*4)
    vk_d = [1e-6]*(4096*4)
    for l in range(d*d):
        if (l%d) >= 0 and (l%d) < (d-k): vk[l] = 1
        if (d-k) <= (l%d) and (l%d) < d: vk_d[l] = 1

    return ((ct<<k)*vk) + ((ct<<(d*d+k-d))*vk_d)

def get_bk(ct, k, d):
    return ct<<(d*k)

def efficient_matmul(mat1, mat2, d):
    # Step 1: get transformed version of mat1 & mat2
    cta_zero = get_cta(mat1, d)
    ctb_zero = get_ctb(mat2, d)
    cta_zero += cta_zero >> d*d
    ctb_zero += ctb_zero >> d*d

    # Step 2: 
    cta_k = []
    ctb_k = []
    for k in range(1, d):
        cta_k.append(get_ak(cta_zero, k, d))
        ctb_k.append(get_bk(ctb_zero, k, d))

    # Step 3:
    ct_ab = cta_zero*ctb_zero # Algorithm 2 - Line 7
    for k in range(0, d-1):
        ct_ab += (cta_k[k]*ctb_k[k])
    return ct_ab


def get_next_length(mat1, mat2, d):

    identity_mat = [0.0]*(4096*4)
    for i in range(d):
        identity_mat[i*d+i] = 1.0
    
    return efficient_matmul(mat1+identity_mat, mat2, d)

def graphanalyticmatmul(graph):
    global N
    
    filt = [1.0]*(N*N)
    while len(filt) < 4096*4:
        filt.append(0.0)
    filt = py_to_eva(filt)
    adjacency_k = graph*filt
    adjacency = (graph*(filt>>N*N))<<(N*N)
    adjacency_k += (adjacency_k>>N*N)  # Matrix is copied for rotations in the algorithm taken from paper
    adjacency += (adjacency>>N*N)  # Matrix is copied for rotations in the algorithm taken from paper

    a_k = get_next_length(adjacency_k, adjacency, N)
    # Commented because of bit modulus error
    # for i in range(3, N):
    #     a_k = get_next_length(a_k, graph, N)

    return a_k

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

def graphanaliticnaivematmul(graph):
    reval = graph
    while len(reval) < 4096*4:
        reval.append(0.0)
    identity = []
    for i in range(4):
        for j in range(4):
            if i==j: identity.append(1.0)
            else: identity.append(0.0)
    while len(identity) < 4096*4:
        identity.append(0.0)

    global N
    n=N
    for row in range(1):
        for col in range(n):
            temp = row_times_col(reval, graph, row, col)
            mask = []
            for mask_ind in range(4096*4):
                if mask_ind == row*n + col : mask.append(1.0)
                else: mask.append(0.0)
            reval += temp*mask

    return reval

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
def simulate_matmul(n):
    m = 4096*4
    global N 
    N = n
    print("Will start simulation for ", n)
    config = {}
    config['warn_vec_size'] = 'false'
    config['lazy_relinearize'] = 'true'
    config['rescaler'] = 'always'
    config['balance_reductions'] = 'true'
    inputs = prepareInput(n, m, p=0.05, for_matmul=True)
    graphanaltic = EvaProgramDriver("graphanaltic", vec_size=m,n=n)
    
    old_head = inputs["Graph"][:N] # Used for early stop by detecting if there is change in the root node at adj matrix
    adj_mat = inputs["Graph"][:N*N]
    time_dict = {"compile" : [], "keygen" : [], "enc": [], "exec" : [], "dec" : [], "ref" : [], "mse" : []}
    for _ in range(N-1):
        with graphanaltic:
            graph = Input('Graph')
            reval = graphanalyticmatmul(graph)
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
        mse = valuation_mse(outputs, reference)  # since CKKS does approximate computations, this is an important measure that depicts the amount of error
        
        time_dict["compile"].append(compiletime)
        time_dict["keygen"].append(keygenerationtime)
        time_dict["enc"].append(encryptiontime)
        time_dict["exec"].append(executiontime)
        time_dict["dec"].append(decryptiontime)
        time_dict["ref"].append(referenceexecutiontime)
        time_dict["mse"].append(mse)

        ret_val = np.array(outputs["ReturnedValue"])

        ret_val = [1 if ret_val[i]>0.4 else 0 for i in range(ret_val.shape[0])]
        unchanged = 0
        for i in range(N):
            if ret_val[i] > 0.5 and old_head[i] > 0.5: unchanged += 1
            elif ret_val[i] < 0.5 and old_head[i] < 0.5: unchanged +=1
        if unchanged == N:
            break
        else:
            # Concat a_k and adjacency matrix
            old_head = ret_val[:N]
            reprepared_input = ret_val[:N*N]
            reprepared_input = np.concatenate((reprepared_input, np.array(adj_mat[:N*N])))
            reprepared_input = reprepared_input.tolist()
            while len(reprepared_input) < 4096*4: reprepared_input.append(0.0)
            inputs["Graph"] = reprepared_input
    print("Connected", np.count_nonzero(ret_val[:N]), "Non-connected:", N-np.count_nonzero(ret_val[:N]))

    compiletime = np.mean(np.array(time_dict["compile"]))
    keygenerationtime =  np.mean(np.array(time_dict["keygen"]))
    encryptiontime = np.mean(np.array(time_dict["enc"]))
    executiontime = np.mean(np.array(time_dict["exec"]))
    decryptiontime = np.mean(np.array(time_dict["dec"]))
    referenceexecutiontime = np.mean(np.array(time_dict["ref"]))
    mse = np.mean(np.array(time_dict["mse"]))

    compiletimes = np.sum(np.array(time_dict["compile"]))
    keygenerationtimes =  np.sum(np.array(time_dict["keygen"]))
    encryptiontimes = np.sum(np.array(time_dict["enc"]))
    executiontimes = np.sum(np.array(time_dict["exec"]))
    decryptiontimes = np.sum(np.array(time_dict["dec"]))
    referenceexecutiontimes = np.sum(np.array(time_dict["ref"]))
    mses = np.sum(np.array(time_dict["mse"]))

    return compiletime, keygenerationtime, encryptiontime, executiontime, decryptiontime, referenceexecutiontime, mse, compiletimes, keygenerationtimes, encryptiontimes, executiontimes, decryptiontimes, referenceexecutiontimes, mses


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
    inputs = prepareInput(n, m, p=0.05)
    graphanaltic = EvaProgramDriver("graphanaltic", vec_size=m,n=n)
    old_head = inputs["Graph"][:N]
    
    time_dict = {"compile" : [], "keygen" : [], "enc": [], "exec" : [], "dec" : [], "ref" : [], "mse" : []}
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
        mse = valuation_mse(outputs, reference)  # since CKKS does approximate computations, this is an important measure that depicts the amount of error
        
        time_dict["compile"].append(compiletime)
        time_dict["keygen"].append(keygenerationtime)
        time_dict["enc"].append(encryptiontime)
        time_dict["exec"].append(executiontime)
        time_dict["dec"].append(decryptiontime)
        time_dict["ref"].append(referenceexecutiontime)
        time_dict["mse"].append(mse)

        ret_val = np.array(outputs["ReturnedValue"])
        ret_val = [1 if ret_val[i]>0.4 else 0 for i in range(ret_val.shape[0])]
        if(np.allclose(old_head, ret_val[:N])): break
        else: old_head = ret_val[:N] 
        inputs["Graph"] = ret_val
    print("Connected", np.count_nonzero(ret_val[:N]), "Non-connected:", N-np.count_nonzero(ret_val[:N]))

    compiletime = np.mean(np.array(time_dict["compile"]))
    keygenerationtime =  np.mean(np.array(time_dict["keygen"]))
    encryptiontime = np.mean(np.array(time_dict["enc"]))
    executiontime = np.mean(np.array(time_dict["exec"]))
    decryptiontime = np.mean(np.array(time_dict["dec"]))
    referenceexecutiontime = np.mean(np.array(time_dict["ref"]))
    mse = np.mean(np.array(time_dict["mse"]))

    compiletimes = np.sum(np.array(time_dict["compile"]))
    keygenerationtimes =  np.sum(np.array(time_dict["keygen"]))
    encryptiontimes = np.sum(np.array(time_dict["enc"]))
    executiontimes = np.sum(np.array(time_dict["exec"]))
    decryptiontimes = np.sum(np.array(time_dict["dec"]))
    referenceexecutiontimes = np.sum(np.array(time_dict["ref"]))
    mses = np.sum(np.array(time_dict["mse"]))

    return compiletime, keygenerationtime, encryptiontime, executiontime, decryptiontime, referenceexecutiontime, mse, compiletimes, keygenerationtimes, encryptiontimes, executiontimes, decryptiontimes, referenceexecutiontimes, mses


if __name__ == "__main__":
    simcnt = 100  # The number of simulation runs, set it to 3 during development otherwise you will wait for a long time
    # For benchmarking you must set it to a large number, e.g., 100
    # Note that file is opened in append mode, previous results will be kept in the file
    resultfile = open("results.csv", "w")  # Measurement results are collated in this file for you to plot later on
    resultfiles = open("resultss.csv", "w")  # Measurement results are collated in this file for you to plot later on
    resultfile.write("NodeCount,SimCnt,CompileTime,KeyGenerationTime,EncryptionTime,ExecutionTime,DecryptionTime,ReferenceExecutionTime,Mse\n")
    resultfiles.write("NodeCount,SimCnt,CompileTime,KeyGenerationTime,EncryptionTime,ExecutionTime,DecryptionTime,ReferenceExecutionTime,Mse\n")
    resultfile.close()
    resultfiles.close()
    
    print("Simulation campaing started:")
    for nc in range(36,64,4): # Node counts for experimenting various graph sizes
        resultfile = open("results.csv", "a") 
        resultfiles = open("resultss.csv", "a") 
        for i in range(simcnt):
            # Call the simulator
            compiletime, keygenerationtime, encryptiontime, executiontime, decryptiontime, referenceexecutiontime, mse, compiletimes, keygenerationtimes, encryptiontimes, executiontimes, decryptiontimes, referenceexecutiontimes, mses = simulate(nc)
            res = str(nc) + "," + str(i) + "," + str(compiletime) + "," + str(keygenerationtime) + "," +  str(encryptiontime) + "," +  str(executiontime) + "," +  str(decryptiontime) + "," +  str(referenceexecutiontime) + "," +  str(mse) + "\n"
            ress = str(nc) + "," + str(i) + "," + str(compiletimes) + "," + str(keygenerationtimes) + "," +  str(encryptiontimes) + "," +  str(executiontimes) + "," +  str(decryptiontimes) + "," +  str(referenceexecutiontimes) + "," +  str(mses) + "\n"
            resultfile.write(res)
            resultfiles.write(ress)
        resultfile.close()
        resultfiles.close()
