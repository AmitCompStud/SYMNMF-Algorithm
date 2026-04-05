import sys


# euclidean distance between vectors
def distance(vector1, vector2):
    sum = 0
    for i in range(len(vector1)):
        sum += (vector1[i] - vector2[i]) ** 2
    return sum ** 0.5


# returns true if distance between all vectors is less than epsilon
def delta_mu_k(vectors1, vectors2, epsilon=0.0001):
    for i in range(len(vectors1)):
        if distance(vectors1[i], vectors2[i]) >= epsilon:
            return False
    return True


def closest_vector(centroids, vector):
    min_dist = sys.float_info.max
    min_index = None
    for i in range(len(centroids)):
        dist = distance(centroids[i], vector)
        if dist < min_dist:
            min_dist = dist
            min_index = i

    return min_index


def calculate_centroids(vector_lists, first_vector):
    centroids = []
    for i in range(len(vector_lists)):
        if len(vector_lists[i]) == 0:
            centroids.append(first_vector)
        else:
            sum = sum_vectors(vector_lists[i])
            centroids.append([sum[j]/len(vector_lists[i]) for j in range(len(sum))])

    return centroids


def sum_vectors(vectors):
    sum = [0 for i in range(len(vectors[0]))]
    for i in range(len(vectors)):
        for j in range(len(vectors[i])):
            sum[j] += vectors[i][j]
    return sum


def kmeans_clustering(vectors, K, max_iterations):
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


def main():
    if len(sys.argv) == 3:  # given K and Iter
        if not sys.argv[1].isdigit():
            print("Incorrect number of clusters!")
            exit(1)
        K = int(sys.argv[1])
        if not sys.argv[2].isdigit():
            print("Incorrect maximum iteration!")
            exit(1)
        max_iter = int(sys.argv[2])
    elif len(sys.argv) == 2:  # Iter not given so Iter=300
        if not sys.argv[1].isdigit():
            print("Incorrect number of clusters!")
            exit(1)
        K = int(sys.argv[1])
        max_iter = 300
    else:
        print("An Error Has Occurred")
        exit(1)

    line1 = sys.stdin.readline()
    d = line1.count(",") + 1
    vectors = [[float(x) for x in line1.split(",")]]
    for line in sys.stdin.readlines():
        vectors.append([float(x) for x in line.split(",")])

    # checking if K and iter values are valid
    if K <= 1 or K >= len(vectors):
        print("Incorrect number of clusters!")
        exit(1)
    if max_iter <= 1 or max_iter >= 800:
        print("Incorrect maximum iteration!")
        exit(1)

    #Running the algo
    clusters = kmeans_clustering(vectors, K, max_iter)
    s = ""
    for cluster in clusters:
        for i in range(len(cluster)):
            s += f"{cluster[i]:.4f}" + ","
        s = s[:-1] + "\n"
    print(s)
    return 0


if __name__ == '__main__':
    main()