import numpy as np
import scipy.io as si
import random
import pickle
import multiprocessing
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from scipy.sparse import csc_matrix, diags
from deap import base
from deap import creator
from deap import tools

def save_obj (obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
 
def load_obj (name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

def construct_W(X, **kwargs):
    """
    Construct the affinity matrix W through different ways

    Notes
    -----
    if kwargs is null, use the default parameter settings;
    if kwargs is not null, construct the affinity matrix according to parameters in kwargs

    Input
    -----
    X: {numpy array}, shape (n_samples, n_features)
        input data
    kwargs: {dictionary}
        parameters to construct different affinity matrix W:
        y: {numpy array}, shape (n_samples, 1)
            the true label information needed under the 'supervised' neighbor mode
        metric: {string}
            choices for different distance measures
            'euclidean' - use euclidean distance
            'cosine' - use cosine distance (default)
        neighbor_mode: {string}
            indicates how to construct the graph
            'knn' - put an edge between two nodes if and only if they are among the
                    k nearest neighbors of each other (default)
            'supervised' - put an edge between two nodes if they belong to same class
                    and they are among the k nearest neighbors of each other
        weight_mode: {string}
            indicates how to assign weights for each edge in the graph
            'binary' - 0-1 weighting, every edge receives weight of 1 (default)
            'heat_kernel' - if nodes i and j are connected, put weight W_ij = exp(-norm(x_i - x_j)/2t^2)
                            this weight mode can only be used under 'euclidean' metric and you are required
                            to provide the parameter t
            'cosine' - if nodes i and j are connected, put weight cosine(x_i,x_j).
                        this weight mode can only be used under 'cosine' metric
        k: {int}
            choices for the number of neighbors (default k = 5)
        t: {float}
            parameter for the 'heat_kernel' weight_mode
        fisher_score: {boolean}
            indicates whether to build the affinity matrix in a fisher score way, in which W_ij = 1/n_l if yi = yj = l;
            otherwise W_ij = 0 (default fisher_score = false)
        reliefF: {boolean}
            indicates whether to build the affinity matrix in a reliefF way, NH(x) and NM(x,y) denotes a set of
            k nearest points to x with the same class as x, and a different class (the class y), respectively.
            W_ij = 1 if i = j; W_ij = 1/k if x_j \in NH(x_i); W_ij = -1/(c-1)k if x_j \in NM(x_i, y) (default reliefF = false)

    Output
    ------
    W: {sparse matrix}, shape (n_samples, n_samples)
        output affinity matrix W
    """

    # default metric is 'cosine'
    if 'metric' not in kwargs.keys():
        kwargs['metric'] = 'cosine'

    # default neighbor mode is 'knn' and default neighbor size is 5
    if 'neighbor_mode' not in kwargs.keys():
        kwargs['neighbor_mode'] = 'knn'
    if kwargs['neighbor_mode'] == 'knn' and 'k' not in kwargs.keys():
        kwargs['k'] = 5
    if kwargs['neighbor_mode'] == 'supervised' and 'k' not in kwargs.keys():
        kwargs['k'] = 5
    if kwargs['neighbor_mode'] == 'supervised' and 'y' not in kwargs.keys():
        print ('Warning: label is required in the supervised neighborMode!!!')
        exit(0)

    # default weight mode is 'binary', default t in heat kernel mode is 1
    if 'weight_mode' not in kwargs.keys():
        kwargs['weight_mode'] = 'binary'
    if kwargs['weight_mode'] == 'heat_kernel':
        if kwargs['metric'] != 'euclidean':
            kwargs['metric'] = 'euclidean'
        if 't' not in kwargs.keys():
            kwargs['t'] = 1
    elif kwargs['weight_mode'] == 'cosine':
        if kwargs['metric'] != 'cosine':
            kwargs['metric'] = 'cosine'

    # default fisher_score and reliefF mode are 'false'
    if 'fisher_score' not in kwargs.keys():
        kwargs['fisher_score'] = False
    if 'reliefF' not in kwargs.keys():
        kwargs['reliefF'] = False

    n_samples, n_features = np.shape(X)

    # choose 'knn' neighbor mode
    if kwargs['neighbor_mode'] == 'knn':
        k = kwargs['k']
        if kwargs['weight_mode'] == 'binary':
            if kwargs['metric'] == 'euclidean':
                # compute pairwise euclidean distances
                D = pairwise_distances(X)
                D **= 2
                # sort the distance matrix D in ascending order
                dump = np.sort(D, axis=1)
                idx = np.argsort(D, axis=1)
                # choose the k-nearest neighbors for each instance
                idx_new = idx[:, 0:k+1]
                G = np.zeros((n_samples*(k+1), 3))
                G[:, 0] = np.tile(np.arange(n_samples), (k+1, 1)).reshape(-1)
                G[:, 1] = np.ravel(idx_new, order='F')
                G[:, 2] = 1
                # build the sparse affinity matrix W
                W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
                bigger = np.transpose(W) > W
                W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
                return W

            elif kwargs['metric'] == 'cosine':
                # normalize the data first
                X_normalized = np.power(np.sum(X*X, axis=1), 0.5)
                for i in range(n_samples):
                    X[i, :] = X[i, :]/max(1e-12, X_normalized[i])
                # compute pairwise cosine distances
                D_cosine = np.dot(X, np.transpose(X))
                # sort the distance matrix D in descending order
                dump = np.sort(-D_cosine, axis=1)
                idx = np.argsort(-D_cosine, axis=1)
                idx_new = idx[:, 0:k+1]
                G = np.zeros((n_samples*(k+1), 3))
                G[:, 0] = np.tile(np.arange(n_samples), (k+1, 1)).reshape(-1)
                G[:, 1] = np.ravel(idx_new, order='F')
                G[:, 2] = 1
                # build the sparse affinity matrix W
                W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
                bigger = np.transpose(W) > W
                W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
                return W

        elif kwargs['weight_mode'] == 'heat_kernel':
            t = kwargs['t']
            # compute pairwise euclidean distances
            D = pairwise_distances(X)
            D **= 2
            # sort the distance matrix D in ascending order
            dump = np.sort(D, axis=1)
            idx = np.argsort(D, axis=1)
            idx_new = idx[:, 0:k+1]
            dump_new = dump[:, 0:k+1]
            # compute the pairwise heat kernel distances
            dump_heat_kernel = np.exp(-dump_new/(2*t*t))
            G = np.zeros((n_samples*(k+1), 3))
            G[:, 0] = np.tile(np.arange(n_samples), (k+1, 1)).reshape(-1)
            G[:, 1] = np.ravel(idx_new, order='F')
            G[:, 2] = np.ravel(dump_heat_kernel, order='F')
            # build the sparse affinity matrix W
            W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
            bigger = np.transpose(W) > W
            W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
            return W

        elif kwargs['weight_mode'] == 'cosine':
            # normalize the data first
            X_normalized = np.power(np.sum(X*X, axis=1), 0.5)
            for i in range(n_samples):
                    X[i, :] = X[i, :]/max(1e-12, X_normalized[i])
            # compute pairwise cosine distances
            D_cosine = np.dot(X, np.transpose(X))
            # sort the distance matrix D in ascending order
            dump = np.sort(-D_cosine, axis=1)
            idx = np.argsort(-D_cosine, axis=1)
            idx_new = idx[:, 0:k+1]
            dump_new = -dump[:, 0:k+1]
            G = np.zeros((n_samples*(k+1), 3))
            G[:, 0] = np.tile(np.arange(n_samples), (k+1, 1)).reshape(-1)
            G[:, 1] = np.ravel(idx_new, order='F')
            G[:, 2] = np.ravel(dump_new, order='F')
            # build the sparse affinity matrix W
            W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
            bigger = np.transpose(W) > W
            W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
            return W

    # choose supervised neighborMode
    elif kwargs['neighbor_mode'] == 'supervised':
        k = kwargs['k']
        # get true labels and the number of classes
        y = kwargs['y']
        label = np.unique(y)
        n_classes = np.unique(y).size
        # construct the weight matrix W in a fisherScore way, W_ij = 1/n_l if yi = yj = l, otherwise W_ij = 0
        if kwargs['fisher_score'] is True:
            W = lil_matrix((n_samples, n_samples))
            for i in range(n_classes):
                class_idx = (y == label[i])
                class_idx_all = (class_idx[:, np.newaxis] & class_idx[np.newaxis, :])
                W[class_idx_all] = 1.0/np.sum(np.sum(class_idx))
            return W

        # construct the weight matrix W in a reliefF way, NH(x) and NM(x,y) denotes a set of k nearest
        # points to x with the same class as x, a different class (the class y), respectively. W_ij = 1 if i = j;
        # W_ij = 1/k if x_j \in NH(x_i); W_ij = -1/(c-1)k if x_j \in NM(x_i, y)
        if kwargs['reliefF'] is True:
            # when xj in NH(xi)
            G = np.zeros((n_samples*(k+1), 3))
            id_now = 0
            for i in range(n_classes):
                class_idx = np.column_stack(np.where(y == label[i]))[:, 0]
                D = pairwise_distances(X[class_idx, :])
                D **= 2
                idx = np.argsort(D, axis=1)
                idx_new = idx[:, 0:k+1]
                n_smp_class = (class_idx[idx_new[:]]).size
                if len(class_idx) <= k:
                    k = len(class_idx) - 1
                G[id_now:n_smp_class+id_now, 0] = np.tile(class_idx, (k+1, 1)).reshape(-1)
                G[id_now:n_smp_class+id_now, 1] = np.ravel(class_idx[idx_new[:]], order='F')
                G[id_now:n_smp_class+id_now, 2] = 1.0/k
                id_now += n_smp_class
            W1 = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
            # when i = j, W_ij = 1
            for i in range(n_samples):
                W1[i, i] = 1
            # when x_j in NM(x_i, y)
            G = np.zeros((n_samples*k*(n_classes - 1), 3))
            id_now = 0
            for i in range(n_classes):
                class_idx1 = np.column_stack(np.where(y == label[i]))[:, 0]
                X1 = X[class_idx1, :]
                for j in range(n_classes):
                    if label[j] != label[i]:
                        class_idx2 = np.column_stack(np.where(y == label[j]))[:, 0]
                        X2 = X[class_idx2, :]
                        D = pairwise_distances(X1, X2)
                        idx = np.argsort(D, axis=1)
                        idx_new = idx[:, 0:k]
                        n_smp_class = len(class_idx1)*k
                        G[id_now:n_smp_class+id_now, 0] = np.tile(class_idx1, (k, 1)).reshape(-1)
                        G[id_now:n_smp_class+id_now, 1] = np.ravel(class_idx2[idx_new[:]], order='F')
                        G[id_now:n_smp_class+id_now, 2] = -1.0/((n_classes-1)*k)
                        id_now += n_smp_class
            W2 = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
            bigger = np.transpose(W2) > W2
            W2 = W2 - W2.multiply(bigger) + np.transpose(W2).multiply(bigger)
            W = W1 + W2
            return W

        if kwargs['weight_mode'] == 'binary':
            if kwargs['metric'] == 'euclidean':
                G = np.zeros((n_samples*(k+1), 3))
                id_now = 0
                for i in range(n_classes):
                    class_idx = np.column_stack(np.where(y == label[i]))[:, 0]
                    # compute pairwise euclidean distances for instances in class i
                    D = pairwise_distances(X[class_idx, :])
                    D **= 2
                    # sort the distance matrix D in ascending order for instances in class i
                    idx = np.argsort(D, axis=1)
                    idx_new = idx[:, 0:k+1]
                    n_smp_class = len(class_idx)*(k+1)
                    G[id_now:n_smp_class+id_now, 0] = np.tile(class_idx, (k+1, 1)).reshape(-1)
                    G[id_now:n_smp_class+id_now, 1] = np.ravel(class_idx[idx_new[:]], order='F')
                    G[id_now:n_smp_class+id_now, 2] = 1
                    id_now += n_smp_class
                # build the sparse affinity matrix W
                W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
                bigger = np.transpose(W) > W
                W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
                return W

            if kwargs['metric'] == 'cosine':
                # normalize the data first
                X_normalized = np.power(np.sum(X*X, axis=1), 0.5)
                for i in range(n_samples):
                    X[i, :] = X[i, :]/max(1e-12, X_normalized[i])
                G = np.zeros((n_samples*(k+1), 3))
                id_now = 0
                for i in range(n_classes):
                    class_idx = np.column_stack(np.where(y == label[i]))[:, 0]
                    # compute pairwise cosine distances for instances in class i
                    D_cosine = np.dot(X[class_idx, :], np.transpose(X[class_idx, :]))
                    # sort the distance matrix D in descending order for instances in class i
                    idx = np.argsort(-D_cosine, axis=1)
                    idx_new = idx[:, 0:k+1]
                    n_smp_class = len(class_idx)*(k+1)
                    G[id_now:n_smp_class+id_now, 0] = np.tile(class_idx, (k+1, 1)).reshape(-1)
                    G[id_now:n_smp_class+id_now, 1] = np.ravel(class_idx[idx_new[:]], order='F')
                    G[id_now:n_smp_class+id_now, 2] = 1
                    id_now += n_smp_class
                # build the sparse affinity matrix W
                W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
                bigger = np.transpose(W) > W
                W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
                return W

        elif kwargs['weight_mode'] == 'heat_kernel':
            G = np.zeros((n_samples*(k+1), 3))
            id_now = 0
            for i in range(n_classes):
                class_idx = np.column_stack(np.where(y == label[i]))[:, 0]
                # compute pairwise cosine distances for instances in class i
                D = pairwise_distances(X[class_idx, :])
                D **= 2
                # sort the distance matrix D in ascending order for instances in class i
                dump = np.sort(D, axis=1)
                idx = np.argsort(D, axis=1)
                idx_new = idx[:, 0:k+1]
                dump_new = dump[:, 0:k+1]
                t = kwargs['t']
                # compute pairwise heat kernel distances for instances in class i
                dump_heat_kernel = np.exp(-dump_new/(2*t*t))
                n_smp_class = len(class_idx)*(k+1)
                G[id_now:n_smp_class+id_now, 0] = np.tile(class_idx, (k+1, 1)).reshape(-1)
                G[id_now:n_smp_class+id_now, 1] = np.ravel(class_idx[idx_new[:]], order='F')
                G[id_now:n_smp_class+id_now, 2] = np.ravel(dump_heat_kernel, order='F')
                id_now += n_smp_class
            # build the sparse affinity matrix W
            W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
            bigger = np.transpose(W) > W
            W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
            return W

        elif kwargs['weight_mode'] == 'cosine':
            # normalize the data first
            X_normalized = np.power(np.sum(X*X, axis=1), 0.5)
            for i in range(n_samples):
                X[i, :] = X[i, :]/max(1e-12, X_normalized[i])
            G = np.zeros((n_samples*(k+1), 3))
            id_now = 0
            for i in range(n_classes):
                class_idx = np.column_stack(np.where(y == label[i]))[:, 0]
                # compute pairwise cosine distances for instances in class i
                D_cosine = np.dot(X[class_idx, :], np.transpose(X[class_idx, :]))
                # sort the distance matrix D in descending order for instances in class i
                dump = np.sort(-D_cosine, axis=1)
                idx = np.argsort(-D_cosine, axis=1)
                idx_new = idx[:, 0:k+1]
                dump_new = -dump[:, 0:k+1]
                n_smp_class = len(class_idx)*(k+1)
                G[id_now:n_smp_class+id_now, 0] = np.tile(class_idx, (k+1, 1)).reshape(-1)
                G[id_now:n_smp_class+id_now, 1] = np.ravel(class_idx[idx_new[:]], order='F')
                G[id_now:n_smp_class+id_now, 2] = np.ravel(dump_new, order='F')
                id_now += n_smp_class
            # build the sparse affinity matrix W
            W = csc_matrix((G[:, 2], (G[:, 0], G[:, 1])), shape=(n_samples, n_samples))
            bigger = np.transpose(W) > W
            W = W - W.multiply(bigger) + np.transpose(W).multiply(bigger)
            return W

def lap_score(X, **kwargs):
    """
    This function implements the laplacian score feature selection, steps are as follows:
    1. Construct the affinity matrix W if it is not specified
    2. For the r-th feature, we define fr = X(:,r), D = diag(W*ones), ones = [1,...,1]', L = D - W
    3. Let fr_hat = fr - (fr'*D*ones)*ones/(ones'*D*ones)
    4. Laplacian score for the r-th feature is score = (fr_hat'*L*fr_hat)/(fr_hat'*D*fr_hat)

    Input
    -----
    X: {numpy array}, shape (n_samples, n_features)
        input data
    kwargs: {dictionary}
        W: {sparse matrix}, shape (n_samples, n_samples)
            input affinity matrix

    Output
    ------
    score: {numpy array}, shape (n_features,)
        laplacian score for each feature

    Reference
    ---------
    He, Xiaofei et al. "Laplacian Score for Feature Selection." NIPS 2005.
    """

    # if 'W' is not specified, use the default W
    if 'W' not in kwargs.keys():
        W = construct_W(X)
    else:
        W = kwargs['W']
    # build the diagonal D matrix from affinity matrix W
    D = np.array(W.sum(axis=1))
    L = W
    tmp = np.dot(np.transpose(D), X)
    D = diags(np.transpose(D), [0])
    Xt = np.transpose(X)
    t1 = np.transpose(np.dot(Xt, D.todense()))
    t2 = np.transpose(np.dot(Xt, L.todense()))
    # compute the numerator of Lr
    D_prime = np.sum(np.multiply(t1, X), 0) - np.multiply(tmp, tmp)/D.sum()
    # compute the denominator of Lr
    L_prime = np.sum(np.multiply(t2, X), 0) - np.multiply(tmp, tmp)/D.sum()
    # avoid the denominator of Lr to be 0
    D_prime[D_prime < 1e-12] = 10000

    # compute laplacian score for all features
    score = 1 - np.array(np.multiply(L_prime, 1/D_prime))[0, :]
    return np.transpose(score)

CXPB = 0.7 # Crossover Rate: the probability to apply Crossover on two selected individuals
POP_SIZE = 1000 # The number of individuals
EPOCHS = 500
FITNESSTYPE = 1 # 1 - Fitness = RF Accuracy | 3 - Fitness = RF Accuracy + Pathway-Pathway similarity + Laplacian score
# a, b and c define the linear combination of accuracy, pathway-pathway similarity and Laplacian score for our fitness function
a = 0.6 # Weight of accuracy
b = 0.2 # Weight of pathway-pathway Jaccard similarity
c = 0.2 # Weight of Laplacian score of selected pathways

dataset = np.load("dataset.npy")
labels = np.load("labels.npy")
jaccardMatrix = np.load("jaccard.npy") # Pathway-pathway similarity based on Jaccard index (ranging from 0 to 1)

# We find maximum similarity value in our pathway-pathway Jaccard matrix to scale our fitness values between 0 and 1
# Maximum similarity is equal to similarity value considering all pathways (i.e. considering all features)
maxJaccardSimilarity = 0
for t in range(0, dataset.shape[1]):
    for v in range(1, dataset.shape[1]):
        if t < v: # Jaccard matrix is specular, so Matrix[i,j] = Matrix[j,i]. We don't consider each value twice!
            maxJaccardSimilarity += jaccardMatrix[t, v]

# We split our dataset with a 80/20 split: 80% for training set and remaining 20% for test set.
# We further split our training set in two equal parts: first one for feature selection with GA; second one for classifiers tuning.

# First Step: 80/20 split of the data, with 20% left for test and 80% for training
train, test, trainLabels, testLabels = train_test_split(dataset, labels, test_size=0.2, random_state=42, stratify = labels)
# Second Step: Splitting training data in two equal part: first one to make individuals evolve and choose the best one (i.e. select best features) and second one for classifiers tuning
trainFeature, trainParams, trainFeatureLabels, trainParamsLabels = train_test_split(train, trainLabels, test_size=0.5, random_state=42, stratify = trainLabels)

# We create our Fitness class. We have to maximize a single objective fitness, that's why we have a single value for weights. 
# Our objective fitness is Random Forest accuracy over 40% of the data
creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax) # We create our class representing an Individual

toolbox = base.Toolbox()

# The function to generate chromosome values (if individual[i] = 1, consider i-th feature; do not if individual[i] = 0)
toolbox.register("attr_bool", random.randint, 0, 1)
# Here we specify how to create an individual. dataset.shape[1] is the number of chromosomes (the number of pathways, i.e. number of features)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, dataset.shape[1])
# Population is an ensemble of individuals!
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# Our function to evaluate fitness value of an individual
def evalAccuracy(individual):
    skf = StratifiedKFold(n_splits=3, shuffle=True)
    indices = [i for i, x in enumerate(individual) if x == 1]
    scoreCV = 0

    for trainCV, testCV in skf.split(trainFeature, trainFeatureLabels):
    	clf = RandomForestClassifier(n_estimators = 301, min_samples_leaf = 1, bootstrap=True, oob_score = True, n_jobs = -1)
        
        trainAfterSelection = trainFeature[trainCV]
        trainAfterSelection = trainAfterSelection[:,indices] # Selecting features according to the chromosomes of the individual
        testAfterSelection = trainFeature[testCV]
        testAfterSelection = testAfterSelection[:, indices]

        # Computing Laplacian Score on train folds only!
        laplacianScores = lap_score(trainFeature[trainCV])
        maxLaplacianScore = sum(laplacianScores) # Maximum laplacian score is the sum of the scores of all our features (i.e. when we do no feature selection)

        # Calculating pathway-pathway similarity and the laplacian score of the solution identified by the current chromosome
        pathSimilarity = 0
        laplacianTotalScore = 0
        for t in range(0, len(indices)):
            laplacianTotalScore += laplacianScores[indices[t]]
            for v in range(1, len(indices)):
                if t < v: # Jaccard matrix is specular, so Matrix[i,j] = Matrix[j,i]. We don't consider each value twice!
                    pathSimilarity += jaccardMatrix[indices[t], indices[v]]
        pathSimilarity /= maxJaccardSimilarity
        laplacianTotalScore /= maxLaplacianScore

        clf.fit(trainAfterSelection, trainFeatureLabels[trainCV])
        scoreCV += accuracy_score(trainFeatureLabels[testCV], clf.predict(testAfterSelection))
    scoreCV /= 3
    if FITNESSTYPE == 1:
        fitnessScore = scoreCV
    elif FITNESSTYPE == 3:
        fitnessScore = a * scoreCV + b * pathSimilarity + c * laplacianTotalScore # Fitness = 3-Fold CV Accuracy, Pathway Similarity and Laplacian Score
    return [fitnessScore]

def selElitistAndTournament(individuals, k_elitist, k_tournament, tournsizeTour):
    return tools.selBest(individuals, k_elitist) + tools.selTournament(individuals, k_tournament, tournsize=tournsizeTour)

toolbox.register("evaluate", evalAccuracy)
toolbox.register("mate", tools.cxOnePoint) # Single point crossover function
# Mutation function. indpb is the indipendent probability to change the value of a chromosome
toolbox.register("mutate", tools.mutFlipBit, indpb=0.3)
# 10% of individuals with the best fitness values were selected to be passed to the next generation.
# After we practise elitism on 10% of our population, we select remaining 90% through tournament technique
toolbox.register("select", selElitistAndTournament, k_elitist=int(0.1*POP_SIZE), k_tournament=POP_SIZE - int(0.1*POP_SIZE), tournsizeTour=int(0.1*POP_SIZE))

pool = multiprocessing.Pool()
toolbox.register("map", pool.map) # Replacing map function with a parallel map impementation

pop = toolbox.population(n=POP_SIZE)
fitnesses = list(toolbox.map(toolbox.evaluate, pop))

for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit
fits = [ind.fitness.values[0] for ind in pop]
g = 0

while max(fits) < 1.0 and g < EPOCHS:
    g = g + 1
    print("-- Generation %i --" % g)
    offspring = toolbox.select(pop)
    # We have to copy our population before we apply crossover and/or mutation.
    # The toolbox.clone() method ensure that we don’t use a reference to the individuals but an completely independent 
    # instance. This is of utter importance since the genetic operators in toolbox will modify the provided objects 
    # in-place
    offspring = list(toolbox.map(toolbox.clone, offspring))
    
    # Crossover
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXPB:
            toolbox.mate(child1, child2)
            # After we applied Crossover, we have a new fitness value for our two individuals
            # We delete old values, in order to find them invalid later on
            del child1.fitness.values
            del child2.fitness.values

    # Mutation
    for mutant in offspring:
        toolbox.mutate(mutant)
        del mutant.fitness.values
    
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind) # Recalculate fitness of invalid (evolved) individuals
    
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    pop[:] = offspring # Saving new generation

    fits = [ind.fitness.values[0] for ind in pop] # Fitness values of new generation
        
    length = len(pop)
    mean = sum(fits) / length
    sum2 = sum(x*x for x in fits)
    std = abs(sum2 / length - mean**2)**0.5
    
    print("  Min %s" % min(fits))
    print("  Max %s" % max(fits))
    print("  Avg %s" % mean)
    print("  Std %s" % std)

# Printing best solution
best_ind = tools.selBest(pop, 1)[0] # Picking best individual
print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))

# Printing results of best individual
indices = [i for i, x in enumerate(best_ind) if x == 1] # Taking indices of pathways selected through Genetic Algorithm in the best individual
skf = StratifiedKFold(n_splits=3, shuffle=True)
scoreCV = 0

# Saving train data both with all features and with selected features only. We do this to plot data with MDS to know their separability
si.savemat('trainNoFeatSelection.mat', dict(data=train, labelsData=trainLabels)) # Train data before feature selection
if FITNESSTYPE == 1:
    si.savemat('trainGA-ACC.mat', dict(data=train[:, indices], labelsData=trainLabels)) # Train data after feature selection
    save_obj(best_ind, 'bestIndividual-ACC')
elif FITNESSTYPE == 3:
    si.savemat('trainGA-ACC_PAT_LAP-' + str(a * 100) + '_' + str(b * 100) + '_' + str(c * 100) + '.mat', dict(data=train[:, indices], labelsData=trainLabels)) # Train data after feature selection
    save_obj(best_ind, 'bestIndividual-' + str(a * 100) + '_' + str(b * 100) + '_' + str(c * 100))

for trainCV, testCV in skf.split(trainFeature, trainFeatureLabels):
    clf = RandomForestClassifier(n_estimators = 301, min_samples_leaf = 1, bootstrap=True, oob_score = True, n_jobs = -1)
    
    trainAfterSelection = trainFeature[trainCV]
    trainAfterSelection = trainAfterSelection[:,indices] # Selecting features according to the chromosomes of the individual
    testAfterSelection = trainFeature[testCV]
    testAfterSelection = testAfterSelection[:, indices]

    # Computing Laplacian Score on train folds only!
    laplacianScores = lap_score(trainFeature[trainCV])
    maxLaplacianScore = sum(laplacianScores) # Maximum laplacian score is the sum of the scores of all our features (i.e. when we do no feature selection)

    # Calculating pathway-pathway similarity and the laplacian score of the solution identified by the current chromosome
    pathSimilarity = 0
    laplacianTotalScore = 0
    for t in range(0, len(indices)):
        laplacianTotalScore += laplacianScores[indices[t]]
        for v in range(1, len(indices)):
            if t < v: # Jaccard matrix is specular, so Matrix[i,j] = Matrix[j,i]. We don't consider each value twice!
                pathSimilarity += jaccardMatrix[indices[t], indices[v]]
    pathSimilarity /= maxJaccardSimilarity
    laplacianTotalScore /= maxLaplacianScore

    clf.fit(trainAfterSelection, trainFeatureLabels[trainCV])
    scoreCV += accuracy_score(trainFeatureLabels[testCV], clf.predict(testAfterSelection))
scoreCV /= 3

print("- - - - Best individual - - - -")
print("Accuracy: %s, Pathway Similarity: %s, Laplacian Score: %s" % (str(scoreCV), str(pathSimilarity), str(laplacianTotalScore)))


# CLASSIFIERS TUNING
trainAfterSelection = trainParams[:, indices]
svm = SVC(C=1.0, kernel='linear')

# SVM

C = [0.001, 0.01, 0.1, 1, 10, 100, 1000] # Choices of SVM C parameter

param_grid = dict(C = C)

estimator = GridSearchCV(svm, cv=3, param_grid=param_grid, n_jobs = -1, verbose = 1)
estimator.fit(trainAfterSelection, trainParamsLabels)

print("The best parameters are %s with a score of %0.4f" % (estimator.best_params_, estimator.best_score_))

# RANDOM FOREST

# Leaf size decision

iterations = 1000
sizes = [1, 2, 5, 7]
averages = [0 for i in range(len(sizes))]

for j in range(iterations):
    print(j)
    for i in range(len(sizes)):
        rf = RandomForestClassifier(bootstrap=True, min_samples_leaf=sizes[i], n_estimators=351, n_jobs=-1, oob_score=True)
        rf= rf.fit(trainAfterSelection, trainParamsLabels)
        averages[i] += rf.oob_score_;
        
averagesFinal = [x/iterations for x in averages]
print(averagesFinal)

# Trees number decision

iterations = 1000
sizes = [51, 101, 151, 251, 301, 351, 501]
averages = [0 for i in range(len(sizes))]

for j in range(iterations):
    print(j)
    for i in range(len(sizes)):
        rf = RandomForestClassifier(bootstrap=True, min_samples_leaf=1, n_estimators=sizes[i], n_jobs=-1, oob_score=True)
        rf= rf.fit(trainAfterSelection, trainParamsLabels)
        averages[i] += rf.oob_score_;

averagesFinal = [x/iterations for x in averages]
print(averagesFinal)

# Random Forest Validation

iterations = 1000
accuracy = []

for i in range(iterations):
    print(i)
    classifier = RandomForestClassifier(bootstrap=True, min_samples_leaf=1, n_estimators=501, n_jobs=-1, oob_score=True)
    classifier.fit(trainAfterSelection, trainParamsLabels)
    accuracy.append(classifier.oob_score_)
print(np.mean(accuracy))

print(accuracy_score(testLabelsFinal, classifier.predict(testFinal)))

# TESTING

svm = SVC(C=1.0)
rf = RandomForestClassifier(bootstrap=True, min_samples_leaf=1, n_estimators=501, n_jobs=-1)

trainAfterSelection = train[:, indices]
testAfterSelection = test[:, indices]

svm.fit(trainAfterSelection, trainLabels)
rf.fit(trainAfterSelection, trainLabels)

print("SVM ACCURACY")
print(accuracy_score(testLabels, svm.predict(testAfterSelection)))

print("RANDOM FOREST ACCURACY")
print(accuracy_score(testLabels, rf.predict(testAfterSelection)))