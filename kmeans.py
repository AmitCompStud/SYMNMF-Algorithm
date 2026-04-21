import sys


def distance(vector1, vector2):
    """
    Calculates the Euclidean distance between two vectors.
    Input:
    vector1: First vector (list of floats)
    vector2: Second vector (list of floats)
    Returns euclidean distance (float)
    """
    sum = 0
    for i in range(len(vector1)):
        sum += (vector1[i] - vector2[i]) ** 2
    return sum ** 0.5


def delta_mu_k(vectors1, vectors2, epsilon=0.0001):
    """
    Checks if the distance between all corresponding vectors is less than epsilon.
    Input:
    vectors1: List of vectors
    vectors2: List of vectors (same length as vectors1)
    epsilon: Threshold for convergence (default 0.0001)
    Returns True if all distances < epsilon, False otherwise
    """
    for i in range(len(vectors1)):
        if distance(vectors1[i], vectors2[i]) >= epsilon:
            return False
    return True


def closest_vector(centroids, vector):
    """
    Finds the index of the closest centroid to a given vector.
    Input:
    centroids: List of centroid vectors
    vector: Target vector
    Returns index of the closest centroid
    """
    min_dist = sys.float_info.max
    min_index = None
    for i in range(len(centroids)):
        dist = distance(centroids[i], vector)
        if dist < min_dist:
            min_dist = dist
            min_index = i

    return min_index


def calculate_centroids(vector_lists, first_vector):
    """
    Calculates new centroids from clusters of vectors.
    Input:
    vector_lists: List of lists, where each sublist contains vectors in a cluster
    first_vector: Fallback vector for empty clusters
    Returns list of new centroid vectors
    """
    centroids = []
    for i in range(len(vector_lists)):
        if len(vector_lists[i]) == 0:
            centroids.append(first_vector)
        else:
            sum = sum_vectors(vector_lists[i])
            centroids.append([sum[j]/len(vector_lists[i]) for j in range(len(sum))])

    return centroids


def sum_vectors(vectors):
    """
    Sums a list of vectors element-wise.
    Input:
    vectors: List of vectors (all same length)
    Returns sum vector (list of floats)
    """
    sum = [0 for i in range(len(vectors[0]))]
    for i in range(len(vectors)):
        for j in range(len(vectors[i])):
            sum[j] += vectors[i][j]
    return sum


def kmeans_clustering(vectors, K, max_iterations):
    """
    Performs K-means clustering on a set of vectors.
    Input:
    vectors: List of data vectors
    K: Number of clusters
    max_iterations: Maximum number of iterations
    Returns final centroids (list of vectors)
    """
    centroids = [vectors[i] for i in range(K)]
    vector_lists = [[] for _ in range(K)]
    for vector in vectors:
        index = closest_vector(centroids, vector)
        vector_lists[index].append(vector)
    new_centroids = calculate_centroids(vector_lists, vectors[0])
    counter = 1
    while not delta_mu_k(new_centroids, centroids) and counter < max_iterations:
        centroids = new_centroids
        #reseting the vector lists
        vector_lists = [[] for _ in range(K)]
        for vector in vectors:
            index = closest_vector(centroids, vector)
            vector_lists[index].append(vector)
        new_centroids = calculate_centroids(vector_lists, vectors[0])
        counter += 1
    return new_centroids
