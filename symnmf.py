import sys
import numpy as np
import symnmfmodule


def error():
    print("An Error Has Occurred")
    sys.exit(1)


def print_matrix(matrix):
    rows = len(matrix)
    for i in range(rows):
        formatted_row = [f"{val:.4f}" for val in matrix[i]]
        print(",".join(formatted_row))


def initialize_H(W, N, k):
    W_average = np.mean(W)
    upper = 2 * np.sqrt(W_average / k)
    H = np.random.uniform(0, upper, size=(N,k))
    return H.tolist()

def handle_goal(goal, matrix, k):
    N = len(matrix)
    if goal == "sym":
        res = symnmfmodule.sym(matrix)
    elif goal == "ddg":
        res = symnmfmodule.ddg(matrix)
    elif goal == "norm":
        res = symnmfmodule.norm(matrix)
    elif goal == "symnmf":
        W = symnmfmodule.norm(matrix) # creating W
        H = initialize_H(W, N, k)
        res = symnmfmodule.symnmf(W, H)
    else:
        error()
    
    return res


def main():
    #setting the seed
    np.random.seed(1234)

    # checking if argv is corect
    if len(sys.argv) != 4:
        error()
    if not sys.argv[1].isdigit(): # if k is not a number, give an error
        error()
    k = int(sys.argv[1])
    goal = sys.argv[2]
    if goal not in ["symnmf", "sym", "ddg", "norm"]: # if goal is not one of the legal options
        error()
    # reading the input file
    file_path = sys.argv[3]
    try:
        input_matrix = np.loadtxt(file_path, delimiter=",")
    except OSError:
        error()
    input_matrix = input_matrix.tolist() # turning into python list
    
    # checking that k is legal
    if k >= len(input_matrix):
        error()
        
    # handle goal
    res = handle_goal(goal, input_matrix, k)
    print_matrix(res)
    return 0
    

if __name__ == "__main__":
    
    main()