import kmeans
import symnmf
import sys
from sklearn.metrics import silhouette_score
import numpy as np


def assign_kmeans(data_points, centroids):
    kmeans_assignments = []
    for point in data_points:
        # getting the distances from current point to all centroids
        distances = [np.linalg.norm(point - c) for c in centroids]
        
        # finding closest cluster index
        closest_cluster = np.argmin(distances)
        kmeans_assignments.append(closest_cluster)
    return kmeans_assignments


def main():
    if len(sys.argv) != 3: # wrong amount of inputs
        symnmf.error()
    
    if not sys.argv[1].isdigit(): # if k is not a number, give an error
        symnmf.error()
    k = int(sys.argv[1])

    file_path = sys.argv[2]
    try:
        data_points = np.loadtxt(file_path, delimiter=",")
    except OSError:
        symnmf.error()
    input_matrix = data_points.tolist() # turning into python list
    #getting the symnmf and kmeans result
    H = symnmf.handle_goal("symnmf", input_matrix, k) # running symnmf algo

    centroids = kmeans.kmeans_clustering(input_matrix, k, 300) # running kmeans algo

    H_np = np.array(H)
    symnmf_assignments = np.argmax(H_np, axis=1) # assigns index of cluster for each data point

    centroids_np = np.array(centroids)
    kmeans_assignments = assign_kmeans(data_points, centroids_np)

    # silhouette_score for each algorithm
    score_kmeans = silhouette_score(data_points, kmeans_assignments)
    score_symnmf = silhouette_score(data_points, symnmf_assignments)
    print(f"nmf: {score_symnmf:.4f}\nkmeans: {score_kmeans:.4f}")
    
    return 0

if __name__ == '__main__':
    np.random.seed(1234)
    main()