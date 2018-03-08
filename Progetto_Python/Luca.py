import numpy as np
import random
import multiprocessing
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from deap import base
from deap import creator
from deap import tools

CXPB = 0.7 # Crossover Rate: the probability to apply Crossover on two selected individuals
POP_SIZE = 1000 # The number of individuals
EPOCHS = 500

dataset = np.load("dataset.npy")
labels = np.load("labels.npy")

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
	clf = RandomForestClassifier(n_estimators = 301, min_samples_leaf = 1, bootstrap=True, oob_score = True, n_jobs = -1)
    indices = [i for i, x in enumerate(individual) if x == 1]
    trainAfterSelection = trainFeature[:, indices] # Selecting features according to the chromosomes of the individual
    return cross_val_score(clf, trainAfterSelection, trainFeatureLabels, cv = 3) # Fitness = 3-Fold CV Accuracy

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