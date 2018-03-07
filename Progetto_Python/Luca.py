import numpy as np
import random
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from deap import base
from deap import creator
from deap import tools

CXPB = 0.7
POP_SIZE = 1000

dataset = np.load("dataset.npy")
labels = np.load("labels.npy")

train, test, trainLabels, testLabels = train_test_split(dataset, labels, test_size=0.2, random_state=42, stratify = labels)
trainFeature, trainParams, trainFeatureLabels, trainParamsLabels = train_test_split(train, trainLabels, test_size=0.5, random_state=42, stratify = trainLabels)

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("attr_bool", random.randint, 0, 1)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, dataset.shape[1])
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

def evalAccuracy(individual):
	clf = RandomForestClassifier(n_estimators = 301, min_samples_leaf = 1, bootstrap=True, oob_score = True, n_jobs = -1)
    indices = [i for i, x in enumerate(individual) if x == 1]
    trainAfterSelection = trainFeature[:, indices] 
    return cross_val_score(clf, trainAfterSelection, trainFeatureLabels)

def selElitistAndTournament(individuals, k_elitist, k_tournament, tournsizeTour):
    return tools.selBest(individuals, k_elitist) + tools.selTournament(individuals, k_tournament, tournsize=tournsizeTour)

toolbox.register("evaluate", evalAccuracy)
toolbox.register("mate", tools.cxOnePoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.3)
toolbox.register("select", selElitistAndTournament, k_elitist=int(0.1*POP_SIZE), k_tournament=POP_SIZE - int(0.1*POP_SIZE), tournsizeTour=(0.1*POP_SIZE))

pop = toolbox.population(n=POP_SIZE)
fitnesses = list(map(toolbox.evaluate, pop))

for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit
fits = [ind.fitness.values[0] for ind in pop]
g = 0

while max(fits) < 1.0 and g < 500:
    g = g + 1
    print("-- Generation %i --" % g)
    offspring = toolbox.select(pop, len(pop))
    offspring = list(map(toolbox.clone, offspring))
    
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXPB:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

    for mutant in offspring:
        toolbox.mutate(mutant)
        del mutant.fitness.values
    
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = map(toolbox.evaluate, invalid_ind)
    
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    pop[:] = offspring

    fits = [ind.fitness.values[0] for ind in pop]
        
    length = len(pop)
    mean = sum(fits) / length
    sum2 = sum(x*x for x in fits)
    std = abs(sum2 / length - mean**2)**0.5
    
    print("  Min %s" % min(fits))
    print("  Max %s" % max(fits))
    print("  Avg %s" % mean)
    print("  Std %s" % std)