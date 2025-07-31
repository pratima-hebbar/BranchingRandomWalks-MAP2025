import math
import random
import matplotlib.pyplot as plot
import matplotlib.colors as color
import functools as ft
import igraph as ig
from matplotlib.colors import LinearSegmentedColormap

def clamped(a, b, c):
    """Return if the value a is clamped between values b and c."""
    return (a > b) & (a < c)

def lerp(a, b, t):
    """Return the linear interpolation of floats a and b with parameter t."""
    return a * (1 - t) + b * t

def get_angles(source, target, sRad, tRad, probability, rangeArr):
    """Identify the angles of which circles at source and target with radii sRad and tRad overlap.
       Then, append these values as ranges to rangeArr. If necessary, split the range into two."""
    d = distBetween(source, target)
    l = (sRad * sRad - tRad * tRad + d * d) / (2 * d)
    h = math.sqrt(sRad*sRad-l*l)
    sharedX = l / d * (target[0]-source[0])
    sharedY = l / d * (target[1]-source[1])
    pmX = h / d * (target[1]-source[1])
    pmY = h / d * (target[0]-source[0])
    low = math.atan2(sharedY-pmY,sharedX+pmX)
    high = math.atan2(sharedY+pmY,sharedX-pmX)
    if high < low:
        rangeArr.append([-math.pi, high, probability])
        rangeArr.append([low, math.pi, probability])
    else:
        rangeArr.append([low, high, probability])

def distBetween(a, b):
    """Get the distance between two positions."""
    return get_scale(arr_dif(a, b))

def get_values(): 
    """Asks the user a series of quesions that they must input an answer to
    and returns an array of their answers.  
    """
    questions = ["How much would you like to weigh the random walk model?\nPlease enter a number greater than zero or enter zero if you do not wish to have this model impact the result.", 
                 "How much would you like to weigh the density model?\nPlease enter a number greater than zero or enter zero if you do not wish to have this model impact the result.", 
                 "How much would you like to weigh the cloud model?\nPlease enter a number greater than zero or enter zero if you do not wish to have this model impact the result.", 
                 "How much would you like to weigh the cluster model?\nPlease enter a number greater than zero ot enter zero if you do not wish to have this model impact the result.", 
                 "How much would you like to weigh the lattice model?\nPlease enter a number greater than zero or enter zero if you do not wish to have this model impact the result.", 
                 "How much would you like to weigh the bias model?\nPlease enter a number greater than zero or enter zero if you do not wish to have this model impact the result.", 
                 "How much would you like to weigh the parent influenced bias model?\nPlease enter a number greater than zero or enter zero if you do not wish to have this model impact the result.", 
                 "For the parent influenced bias model only: \nWhen calculating a node's position, how many degrees would you like its angle to vary from its parent's angle?\nThe parent is the node that produced the node and its angle is the direction it moved on the graph.\nPlease enter a number in between 0 and 180.", 
                 "When calculating a node's position, how many degrees would you like a node's angle to vary from another node that shares the same parent as it?\nThis applies to all models so if you do not want any variation, please enter zero.\nOtherwise, please enter a number no greater than 180.", 
                 "For the cluster model only: \nWhen calculating a node's position, how many degrees would you like its angle to vary from the guiding angle?\nThis model picks a direction for all nodes within a calculated distance to travel in give or take the input.\nPlease enter a number in between 0 and 180.", 
                 "For the bias model and the parent influenced bias model only: \nWhen calculating a node's position, how many degrees would you like it to varies from the direction it chooses?\nThese models pick one of three directions and have its nodes travel in that direction give or take the input.\nPlease enter a number in between 0 and 180.",
                 "How many kids would you like each node to have?", 
                 "How many generations would you like this program to run for?",
                 "How many of these generations would you like to be random?"]
    # an array of all questions that the user will answer with their input
    answers = []
    # an array of all the users inputed answers
    index = 0
    # the current index of the questions array
    while index < len(questions):
        print(questions[index])
        # prints the question at the current index
        answer = input()
        # saves the users answer
        index += 1
        try:
            num_answer = float(answer)
            # checks the answer to see if it can become a float
            if num_answer < 0:
                raise ValueError()
            answers.append(num_answer)
            # appends it to the answers array if it can
        except (ValueError, TypeError):
            index -= 1
            # reasks the question if it can't
    return answers
    # returns the array containing the user's answers

def get_answers(arr):
    """Determines, based on user input, whether or not a model will be needed,
    creates a true or false array that indicates so, and appends that array to
    the array of user inputs.

    Init Args:
    arr -- The array of user inputs.
    """
    answers = []
    # an array of true or false values depending on what needs to be run
    for i in range(7):
        if arr[i] <= 0:
            answers.append(False)
            # if the value of the weight is less than or equal to zero, the model need not be run
        else:
            answers.append(True)
            # else the model should be run
    return answers + arr
    # returns the two arrays appended together to be read

answers = get_answers(get_values())
# an array of whether or not a model is needed and their weights, among other things
child_count = int(answers[18])
# the number of children each node will have
random_gen = int(answers[20])
# the generation at which the model finishes completely random modeling
if child_count <= 0:
    child_count = 2
    # in case an invalid child count is given

class Simulation:
    """A class dedicated to running simulations of a branching random walk."""
    def __init__(self, **kwargs):
        self.pos = kwargs.get('pos', [0, 0])
        # the position of the node
        self.par_var = math.radians(kwargs.get('par_var', 15))
        # the variance the kid nodes have from their parent node's angle
        self.sib_var = math.radians(kwargs.get('sib_var', 15))
        # the variance the kid nodes have from their sibling node's angle
        self.clu_var = math.radians(kwargs.get('clu_var', 15))
        # the variance the kid nodes have from the guiding angle
        self.bias_var = math.radians(kwargs.get('bias_var', 15))
        # the variance the kid nodes have from the biased angle
        self.gens = None
        # an double array of all the nodes, sorted by generation
        self.root = None
        self.generation = 0
        # the current generation number
        self.included = [kwargs.get('random_walk', False),
                         kwargs.get('density', False),
                         kwargs.get('cloud', False),
                         kwargs.get('cluster', False), 
                         kwargs.get('lattice', False), 
                         kwargs.get('triangle', False), 
                         kwargs.get('bias', False)]
        # an array of whether or not a model will be used
        self.weights = [kwargs.get('random_walk_weight', 1 if self.included[0] else 0),
                        kwargs.get('density_weight', 1 if self.included[1] else 0),
                        kwargs.get('cloud_weight', 1 if self.included[2] else 0),
                        kwargs.get('clustering_weight', 1 if self.included[3] else 0), 
                        kwargs.get('lattice_weight', 1 if self.included[4] else 0), 
                        kwargs.get('triangle_weight', 1 if self.included[5] else 0), 
                        kwargs.get('bias_weight', 1 if self.included[6] else 0)]
        # an array of the weights of each model
        self.graph = None
        # a graph of the nodes (to be constructed later)
        self.clusters = None
        # a collection of the count of clusters (to be constructed later)
        
    def initialize(self):
        """Reset all local parameters to those of an unrun simulation."""
        self.gens = [[]]
        self.clusters = []
        self.root = Node(self.pos, None, 0, 0, self.par_var, self.sib_var, self.clu_var, self.bias_var, self.gens[0], self.weights, self.included)
        self.generation = 0
        self.graph = None
        
    def run_gen(self):
        """Run a generation of the simulation."""
        next_gen = []
        # create a new array for nodes of the upcoming generation.
        num = len(self.gens[self.generation])
        # get the count of the previous generation's nodes.
        edges = []
        # create an array to store the edges.
        for source_index in range(num):
            for target_index in range(source_index + 1, num):
                if get_scale(arr_dif(self.gens[self.generation][source_index].position, self.gens[self.generation][target_index].position)) < f_t(self.generation):
                    edges.append([source_index, target_index])
        # if the edges between two nodes is less than f(t), include it in the graph.
        g = ig.Graph(num, edges)
        # create the graph.
        components = g.connected_components(mode='weak')
        # get the components in the graph.
        self.clusters.append(len(components))
        # record the number of components in the graph.
        angles = []
        # create an array to store the angles for each component.
        for _ in components:
            angles.append(random.random() * 2 * math.pi)
        # for each component, store a random angle.
        for index, component in enumerate(components):
            for node in component:
                self.gens[self.generation][node].set_guiding_angle(angles[index])
        # for each node in each component, assign its guiding angle to that of its component.
        for node in self.gens[self.generation]:
            node.make_children(next_gen)
        # generate the upcoming generation.
        self.generation += 1
        # increase the generation count.
        self.gens.append(next_gen)
        # store the new generation.
        
    def run_gens(self, num):
        """Run a series of generations"""
        for _ in range(num):
            self.run_gen()
            print("Finished gen", self.generation)
            
    def plot_end(self, env):
        for node in self.gens[-1]:
            node.plot_end(env, (0, 0, 0, 1))
        # plots all the nodes in the final generation
    
    def plot_path(self, env):
        for node in self.gens[0]:
            node.plot_path(env, 0, 1)
        # plots all the nodes and their path with rainbow colors

    def plot_other_path(self, env):
        for node in self.gens[0]:
            node.plot_other_path(env)
        # plots all the nodes and their path with a single color

    def heat_end(self, env):
        x = []
        # an array of x-coordinates
        y = []
        # an array of y-coordinates
        for node in self.gens[-1]:
            x.append(node.position[0])
            # appends all of the final generation's nodes x-coordinates
            y.append(node.position[1])
            # appends all of the final generation's nodes y-coordinates
        my_colormap = LinearSegmentedColormap.from_list("my colormap", ['#000000', '#ff0000', '#ffff00', '#00ff00', '#00ffff', '#0000ff', '#ff00ff'], N = 100)
        # creates a colormap used in the histogram
        env.hist2d(x, y, bins = 75, range = [[-10, 10], [-10, 10]], cmap = my_colormap)
        # constructs a histogram of all the final generation's points

    def heat_path(self, env):
        x = []
        # an array of x-coordinates
        y = []
        # an array of y-coordinates
        for thing in self.gens:
            for node in thing:
                x.append(node.position[0])
                # appends all of the nodes x-coordinates
                y.append(node.position[1])   
                # appends all of the final generation's nodes y-coordinates
        my_colormap = LinearSegmentedColormap.from_list("my colormap", ['#000000', '#ff0000', '#ffff00', '#00ff00', '#00ffff', '#0000ff', '#ff00ff'], N = 100)
        # creates a colormap used in the histogram
        env.hist2d(x, y, bins = 75, range = [[-10, 10], [-10, 10]], cmap = my_colormap)
        # constructs a histogram of all points
        

class Node:
    """A class representing a Node in the branching random walk simulation."""
    def __init__(self, pos, par, ang, gen, par_var, sib_var, clu_var, bias_var, population, w, i):
        self.position = pos
        # the position of this node.
        self.parent = par
        # the parent of this node.
        self.angle = ang
        # the angle this node traveled from its parent.
        self.children = []
        # an array of this node's children.
        self.generation = gen
        # the generation this node is a part of.
        self.guiding = 0
        # the guiding angle of this node.
        self.parent_variance = par_var
        # simulation parameter.
        self.sibling_variance = sib_var
        # simulation parameter.
        self.cluster_variance = clu_var
        # simulation parameter.
        self.bias_variance = bias_var
        # simulation parameter.
        self.population = population
        # this node's cousins.
        self.population.append(self)
        # add this node to the arraay of those of the same generation.
        self.weights = w
        # simulation parameter.
        self.included = i
        # simulation parameter.
        
    def set_guiding_angle(self, angle):
        """Set the guiding angle of this node."""
        self.guiding = angle
    
    def make_children(self, next):
        """Make children for this node according to simulation parameters and past nodes."""
        density = density_func(self, self.population)
        # calculate the density about this node.
        rand_arr = [0, 0]
        density_arr = [0, 0]
        cloud_arr = [0, 0]
        cluster_arr = [0, 0]
        lattice_arr = [0, 0]
        triangle_arr = [0, 0]
        bias_arr = [0, 0]
        # create arrays for each model that will have no effect on the final model unless they are assigned a value.
        # note that the following models 
        if self.included[0] or self.generation < random_gen:
            rand_angle = random.random() * 2 * math.pi
            rand_arr = [math.cos(rand_angle), math.sin(rand_angle)]
        # if this model includes the random walk model, assign meaningful value to its corresponding array
        if self.included[1] and self.generation >= random_gen:
            density_angle = random.random() * 2 * math.pi
            density_arr = to_unit(arr_sum(scale_arr([math.cos(self.angle), math.sin(self.angle)], density), [math.cos(density_angle), math.sin(density_angle)]))
        if self.included[2] and self.generation >= random_gen:
            true_region = [[-math.pi, math.pi, 1]]
            regions = []
            cousins = self.get_child_cousins()
            for degree in range(len(cousins)):
                r = radius(self.generation, degree)
                for cousin in cousins[degree]:
                    if clamped(distBetween(self.position, cousin.position), 1 - r, 1 + r):
                        get_angles(self.position, cousin.position, 1, r, 0.1, regions)
            for region in regions:
                j = 0
                while (j < len(true_region)) and (true_region[j][1] < region[1]):
                    if (region[0] < true_region[j][1]):
                        if region[0] <= true_region[j][0]:
                            true_region[j][2] *= region[2]
                            region[0] = true_region[j][1]
                        else:
                            upper = true_region[j][1]
                            true_region[j][1] = region[0]
                            true_region.insert(j + 1, [region[0], upper, true_region[j][2] * region[2]])
                            region[0] = upper
                            j += 1
                    j += 1
                if (j < len(true_region)) & (region[1] > region[0]):
                    if true_region[j][1] == region[1]:
                        true_region[j][2] *= region[2]
                    elif (region[0] > true_region[j][0]):
                        upper = true_region[j][1]
                        true_region[j][1] = region[1]
                        true_region.insert(j + 1, [region[1], upper, true_region[j][2]])
                        true_region[j][2] *= region[2]
            index = random.choices(range(len(true_region)), map(lambda x: x[2] * (x[1]- x[0]) / math.pi / 2, true_region))[0]
            cloud_angle = lerp(true_region[index][0], true_region[index][1], random.random())
            cloud_arr = [math.cos(cloud_angle), math.sin(cloud_angle)]
        if self.included[3] and self.generation >= random_gen:
            cluster_angle = self.guiding + random.random() * 2 * self.cluster_variance - self.cluster_variance
            cluster_arr = [math.cos(cluster_angle), math.sin(cluster_angle)]
        if self.included[4] and self.generation >= random_gen:
            lattice_angles = [0, 90, 180, 270]
            # an array of all possible lattice angles
            lattice_selected = random.randint(0, 3)
            # a random index for the array of lattice angles
            lattice_angle = math.radians(lattice_angles[lattice_selected])
            # the selected angle for the node to travel
            lattice_arr = [math.cos(lattice_angle), math.sin(lattice_angle)]
            # the coordinates of the new node
        if self.included[5] and self.generation >= random_gen:
            rand_val = random.randint(1, 22)
            # a random value that determines which direction the node goes in
            if rand_val < 8:
                triangle_angle = math.radians(0) + random.random() * 2 * self.bias_variance - self.bias_variance
                # an angle within 10 degrees of 0
            elif rand_val < 15:
                triangle_angle = math.radians(120) + random.random() * 2 * self.bias_variance - self.bias_variance
                # an angle within 10 degrees of 120
            elif rand_val < 22:
                triangle_angle = math.radians(240) + random.random() * 2 * self.bias_variance - self.bias_variance
                # an angle within 10 degrees of 240
            else:
                triangle_angle = random.random() * 2 * math.pi
                # a completely random angle
            triangle_arr = [math.cos(triangle_angle), math.sin(triangle_angle)]
            # the coordinates of the new node
        if self.included[6] and self.generation >= random_gen:
            rand_val = random.randint(1, 22)
            # a random value that determines which direction the node goes in
            if self.generation == 0:
                if rand_val < 8:
                    bias_angle = math.radians(0) + random.random() * 2 * self.bias_variance - self.bias_variance
                    # an angle within 10 degrees of 0
                elif rand_val < 15:
                    bias_angle = math.radians(120) + random.random() * 2 * self.bias_variance - self.bias_variance
                    # an angle within 10 degrees of 120
                elif rand_val < 22:
                    bias_angle = math.radians(240) + random.random() * 2 * self.bias_variance - self.bias_variance
                    # an angle within 10 degrees of 240
                else:
                    bias_angle = random.random() * 2 * math.pi
                    # a completely random angle
            else:
                if rand_val < 22:
                    bias_angle = self.angle + random.random() * 2 * self.parent_variance  - self.parent_variance
                    # an angle within 10 degrees of its parent's angle
                else:
                    bias_angle = random.random() * 2 * math.pi
                    # a completely random angle
            bias_arr = [math.cos(bias_angle), math.sin(bias_angle)]
            # the coordinates of the new node
        if self.generation < random_gen:
            delta = to_unit(arr_sum(scale_arr(rand_arr, 1)))
            # only weighs the random model
        else:
            delta = to_unit(arr_sum(scale_arr(rand_arr, self.weights[0]),
                                    scale_arr(density_arr, self.weights[1]),
                                    scale_arr(cloud_arr, self.weights[2]),
                                    scale_arr(cluster_arr, self.weights[3]), 
                                    scale_arr(lattice_arr, self.weights[4]), 
                                    scale_arr(triangle_arr, self.weights[5]), 
                                    scale_arr(bias_arr, self.weights[6])))
        for _ in range(child_count):
            bearing = math.atan2(delta[1], delta[0]) + random.random() * 2 * self.sibling_variance - self.sibling_variance
            self.children.append(Node(arr_sum(self.position, [math.cos(bearing), math.sin(bearing)]), self, bearing, self.generation + 1, self.parent_variance, self.sibling_variance, self.cluster_variance, self.bias_variance, next, self.weights, self.included))
            
    def get_child_cousins(self):
        cousins = [[]]
        uncommon_ancestor = self
        common_ancestor = self.parent
        i = 1
        while common_ancestor is not None:
            cousins.append([])
            for cousin in common_ancestor.children:
                if cousin is not uncommon_ancestor:
                    cousin.get_child_cousins_helper(cousins[i], i)
            uncommon_ancestor = common_ancestor
            common_ancestor = common_ancestor.parent
            i += 1
        return cousins
    
    def get_child_cousins_helper(self, arr, index):
        if index == 1:
            for child in self.children:
                arr.append(child)
        else:
            for child in self.children:
                child.get_child_cousins_helper(arr, index - 1)
    
    def plot_end(self, env, hue):
        env.plot([self.position[0]], [self.position[1]], 'o', color = hue)

    def plot_path(self, env, hue, depth):
        for i, child in enumerate(self.children):
            new_hue = hue + i / (child_count**depth)
            col = color_alpha(color.hsv_to_rgb((new_hue, 1, 1)), 0.1)
            env.plot([self.position[0], child.position[0]], [self.position[1], child.position[1]], color = col)
            child.plot_path(env, new_hue, depth + 1)
        if len(self.children) == 0:
            self.plot_end(env, color_alpha(color.hsv_to_rgb((hue, 1, 1)), 0.1))

    def plot_other_path(self, env):
        for i, child in enumerate(self.children):
            col = "#00ffff10"
            env.plot([self.position[0], child.position[0]], [self.position[1], child.position[1]], color = col)
            child.plot_other_path(env)
        if len(self.children) == 0:
            self.plot_end(env, "#00ffff10")

def arr_sum(*args):
    output = [0, 0]
    for arg in args:
        output[0] += arg[0]
        output[1] += arg[1]
    return output

def arr_dif(arrA, arrB):
    output = []
    for i in range(len(arrA)):
        output.append(arrA[i] - arrB[i])
    return output

def to_unit(arr):
    scale = get_scale(arr)
    if scale < 0.0001:
        return 0
    output = []
    for i in range(len(arr)):
        output.append(arr[i] / scale)
    return output

def get_scale(arr):
    scale = 0
    for i in arr:
        scale += i * i
    return math.sqrt(scale)

def scale_arr(arr, scale):
    output = []
    for i in arr:
        output.append(i * scale)
    return output

def density_func(s, pop):
    count = 0
    for cousin in pop:
        if get_scale(arr_dif(s.position, cousin.position)) <= 1:
            count += 1 / 4
    return count

def color_alpha(col, alpha):
    return (col[0], col[1], col[2], alpha)

def f_t(t):
    return 0.012 * ((2 * t + 1) / (t + 1))

def radius(gen, degree):
    return 0.1

fig, (ax1, ax2, ax3) = plot.subplots(3, 6, figsize = (20, 10))
plot.subplots_adjust(left=.025, bottom=.05, right=.975, top=.9)
axX = [ax1, ax2, ax3]
sim = Simulation(random_walk = answers[0], density = answers[1], cloud = answers[2], cluster = answers[3], lattice = answers[4], triangle = answers[5], bias = answers[6], random_walk_weight = answers[7], density_weight = answers[8], cloud_weight = answers[9], clustering_weight = answers[10], lattice_weight = answers[11], triangle_weight = answers[12], bias_weight = answers[13], par_var = answers[14], sib_var = answers[15], clu_var = answers[16], bias_var = answers[17])
title_string = "Parameters: "
for val in answers:
    title_string += str(val) + ", "
fig.suptitle(title_string)
for a in range(3):
    for b in range(3):
        sim.initialize()
        sim.run_gens(int(answers[19]))
        sim.plot_path(axX[a][b * 2])
        # sim.plot_end(axX[a][b * 2])
        # sim.heat_end(axX[a][b * 2])
        # sim.heat_path(axX[a][b * 2])
        # sim.plot_other_path(axX[a][b * 2])
        # sim.plot_path(axX[a][b * 2 + 1])
        # sim.plot_end(axX[a][b * 2 + 1])
        sim.heat_end(axX[a][b * 2 + 1])
        # sim.heat_path(axX[a][b * 2 + 1])
        # sim.plot_other_path(axX[a][b * 2 + 1])
        axX[a][b * 2].axis([-int(answers[19]), int(answers[19]), -int(answers[19]), int(answers[19])])
        axX[a][b * 2 + 1].axis([-int(answers[19]), int(answers[19]), -int(answers[19]), int(answers[19])])
        axX[a][b * 2].set_title(f"Sim #{a * 3 + b + 1}: Full plot")
        axX[a][b * 2 + 1].set_title(f"Sim #{a * 3 + b + 1}: 2D Histogram")

plot.show()
# Models to recreate:
# Density repulsion
# Cloud repulsion
# Clustering
# True Random
