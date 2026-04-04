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

        H = np.array(res)
        HT = H.T

        # 1. Multiply H and H^T using the @ operator
        H_multiplied = H @ HT
        
        # 2. Make sure W is also a NumPy array for a clean comparison
        W_np = np.array(W)
        
        # 3. Print them out to eyeball the difference!
        print("=== H * H^T Matrix ===")
        print_matrix(H_multiplied.tolist())
        
        print("\n=== Original W Matrix ===")
        print_matrix(W)
        print('\n')
        difference = np.linalg.norm(W_np - H_multiplied)
        print(f"\nFinal Difference (Frobenius Norm): {difference:.6f}")
    else:
        error()
    
    print_matrix(res)


def main():
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
    handle_goal(goal, input_matrix, k)
    

if __name__ == "__main__":
    np.random.seed(1234)
    main()