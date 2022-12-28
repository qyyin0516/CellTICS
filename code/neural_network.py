import numpy as np
import math
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()


def create_placeholders(n_x, n_y):
    X = tf.placeholder(tf.float32, [n_x, None], name="X")
    Y = tf.placeholder(tf.float32, [n_y, None], name="Y")
    return X, Y


def random_mini_batches(X, Y, mini_batch_size=32):
    m = X.shape[1]
    mini_batches = []

    permutation = list(np.random.permutation(m))
    shuffled_X = X[:, permutation]
    shuffled_Y = Y[:, permutation].reshape((Y.shape[0], m))

    num_complete_minibatches = math.floor(m / mini_batch_size)
    for k in range(0, num_complete_minibatches):
        mini_batch_X = shuffled_X[:, k * mini_batch_size: k * mini_batch_size + mini_batch_size]
        mini_batch_Y = shuffled_Y[:, k * mini_batch_size: k * mini_batch_size + mini_batch_size]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    if m % mini_batch_size != 0:
        mini_batch_X = shuffled_X[:, num_complete_minibatches * mini_batch_size: m]
        mini_batch_Y = shuffled_Y[:, num_complete_minibatches * mini_batch_size: m]
        mini_batch = (mini_batch_X, mini_batch_Y)
        mini_batches.append(mini_batch)
    return mini_batches


def initialize_parameters(layers_node, n):
    m = len(layers_node)
    W = {}
    b = {}
    M = {}
    for i in range(n - 1):
        j = m - i - 1
        weight = tf.get_variable("W" + str(i + 1), [len(layers_node[j - 1]), len(layers_node[j])],
                                 initializer=tf.keras.initializers.glorot_normal())
        W["W" + str(i + 1)] = weight
        bias = tf.get_variable("b" + str(i + 1), [len(layers_node[j - 1]), 1], initializer=tf.zeros_initializer())
        b["b" + str(i + 1)] = bias
        mask = tf.get_variable("M" + str(i + 1), [len(layers_node[j - 1]), len(layers_node[j])])
        M["M" + str(i + 1)] = mask
    W["W" + str(n)] = tf.get_variable("W" + str(n), [len(layers_node[0]), len(layers_node[m - n])],
                                 initializer=tf.keras.initializers.glorot_normal())
    b["b" + str(n)] = tf.get_variable("b" + str(n), [len(layers_node[0]), 1], initializer=tf.zeros_initializer())
    parameters = {}
    parameters.update(W)
    parameters.update(b)
    parameters.update(M)
    return parameters


def forward_propagation(X, Y, parameters, masking, n, gamma=0.0001):
    W = {}
    b = {}
    M = {}
    for i in range(n):
        W[i + 1] = parameters["W" + str(i + 1)]
        b[i + 1] = parameters["b" + str(i + 1)]
        if i != (n - 1):
            M[i + 1] = parameters["M" + str(i + 1)]
            M[i + 1] = tf.assign(M[i + 1], masking[len(masking) - i - 1])
    Z = {}
    A = {}
    for i in range(1, n + 1):
        if i == 1:
            Z[i] = tf.add(tf.matmul(tf.multiply(M[i], W[i]), X),  b[i])
            A[i] = tf.nn.tanh(Z[i])
        elif i == n:
            Z[i] = tf.add(tf.matmul(W[i], A[i - 1]), b[i])
            A[i] = Z[i]
        else:
            Z[i] = tf.add(tf.matmul(tf.multiply(M[i], W[i]), A[i - 1]),  b[i])
            A[i] = tf.nn.tanh(Z[i])
    logits = tf.transpose(Z[n])
    labels = tf.transpose(Y)
    reg_term = gamma * tf.nn.l2_loss(W[1])
    for i in range(2, n + 1):
        reg_term = reg_term + gamma * tf.nn.l2_loss(W[i])
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=logits, labels=labels) + reg_term)
    return A, cost


def model(X_train, Y_train, X_test, layers_node, masking, output_layer, learning_rate=0.001, num_epochs=50,
          minibatch_size=32, gamma=0.0001, print_cost=False):
    tf.reset_default_graph()
    (n_x, m) = X_train.shape
    n_y = Y_train.shape[0]
    X, Y = create_placeholders(n_x, n_y)
    parameters = initialize_parameters(layers_node, output_layer)
    A, cost = forward_propagation(X, Y, parameters, masking, output_layer, gamma=gamma)
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    init = tf.global_variables_initializer()

    with tf.Session() as sess:
        sess.run(init)
        for epoch in range(num_epochs):
            epoch_cost = 0.
            num_minibatches = int(m / minibatch_size)
            minibatches = random_mini_batches(X_train, Y_train, minibatch_size)
            for minibatch in minibatches:
                (minibatch_X, minibatch_Y) = minibatch
                _, minibatch_cost = sess.run([optimizer, cost], feed_dict={X: minibatch_X, Y: minibatch_Y})
                epoch_cost += minibatch_cost / num_minibatches
            if print_cost == True:
                print("Cost after epoch %i: %f" % (epoch, epoch_cost))
        output_train = sess.run(A, feed_dict={X: X_train})
        output_test = sess.run(A, feed_dict={X: X_test})
        return output_train, output_test
