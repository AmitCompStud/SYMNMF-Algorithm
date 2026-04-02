import sys
import numpy as np


def error():
    print("An Error Has Occurred")
    sys.exit(1)


def sym():



def main():
    #checking if argv is corect
    if(sys.argv != 4):
        error()
    
    if not sys.argv[1].isdigit(): # if k is not a number, give an error
        error()
    k = int(sys.argv[1])
    goal = argv[2]
    if goal not in ["symnmf", "sym", "ddg", "norm"]: # if goal is not one of the legal options
        error()
    file_path = sys.argv[3]
    
    

        






if __name__ == "__main__":
    main()