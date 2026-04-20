"""Utility classes and functions for the Galacticus normalizing-flow emulator.

Provides coordinate-transformation helpers (:func:`norm_transform` and
:func:`norm_transform_inv`) and a Keras-based Real-valued Non-Volume
Preserving (RealNVP) normalizing flow model (:class:`RealNVP`) used to
emulate the Galacticus semi-analytic model parameter space.
"""
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import regularizers
import numpy as np
import tensorflow_probability as tfp
from scipy import stats

# This function transforms data from normalized coordinates to hypercube coordinates
def norm_transform(data, min_val, max_val):
    """Transform data from normalized coordinates to hypercube coordinates.

    Rescales each column of *data* so that the minimum observed value maps to
    *min_val* and the maximum observed value maps to *max_val*.

    Parameters
    ----------
    data : array_like
        Input data array, shape (n_samples, n_features).
    min_val : float
        Lower bound of the target hypercube interval.
    max_val : float
        Upper bound of the target hypercube interval.

    Returns
    -------
    data_min : ndarray
        Per-feature minimum values of *data*.
    data_max : ndarray
        Per-feature maximum values of *data*.
    transformed : ndarray
        Rescaled data in the range [*min_val*, *max_val*].
    """
    data_min = np.nanmin(data, axis = 0)
    data_max = np.nanmax(data, axis = 0)
    sigma_data = (data - data_min)/(data_max - data_min)
    return data_min, data_max, sigma_data*(max_val - min_val) + min_val

# This function transforms data from hypercube coordinates to normalized coordinates
def norm_transform_inv(norm_data, data_min, data_max, min_val, max_val):
    """Transform data from hypercube coordinates back to normalized coordinates.

    Inverts the rescaling performed by :func:`norm_transform`.

    Parameters
    ----------
    norm_data : array_like
        Data in hypercube coordinates, i.e. values in [*min_val*, *max_val*].
    data_min : array_like
        Per-feature minimum values used when the forward transform was applied.
    data_max : array_like
        Per-feature maximum values used when the forward transform was applied.
    min_val : float
        Lower bound of the hypercube interval.
    max_val : float
        Upper bound of the hypercube interval.

    Returns
    -------
    ndarray
        Data rescaled back to the original (normalized) coordinate range.
    """
    sigma_data = (norm_data - min_val)/(max_val - min_val)
    return sigma_data*(data_max - data_min) + data_min

### CREATING EMULATOR CLASSES/FUNCTIONS ###

# Creating an individual coupling layer with keras API.
def Coupling(input_shape, output_dim = 256, reg = 0.01):
    """Build a single Real-NVP coupling-layer network.

    Constructs a Keras model with two parallel dense sub-networks (the scale
    network *s* and the translation network *t*), each consisting of five
    fully-connected layers with ReLU activations, followed by a ``tanh``
    output layer.  The outputs are used to define an invertible affine
    transformation in the normalizing-flow model.

    Parameters
    ----------
    input_shape : int
        Dimensionality of the input (and output) of the coupling layer.
    output_dim : int, optional
        Width of each hidden dense layer. Defaults to 256.
    reg : float, optional
        L2 regularization coefficient applied to all weight matrices.
        Defaults to 0.01.

    Returns
    -------
    keras.Model
        A Keras model whose outputs are ``[s, t]``, each of shape
        ``(batch, input_shape)``.
    """
    input = keras.layers.Input(shape=(input_shape,))

    t_layer_1 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(input)
    t_layer_2 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(t_layer_1)
    t_layer_3 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(t_layer_2)
    t_layer_4 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(t_layer_3)
    t_layer_5 = keras.layers.Dense(
        input_shape, activation="tanh", kernel_regularizer=regularizers.l2(reg)
    )(t_layer_4)

    s_layer_1 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(input)
    s_layer_2 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(s_layer_1)
    s_layer_3 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(s_layer_2)
    s_layer_4 = keras.layers.Dense(
        output_dim, activation="relu", kernel_regularizer=regularizers.l2(reg)
    )(s_layer_3)
    s_layer_5 = keras.layers.Dense(
        input_shape, activation="tanh", kernel_regularizer=regularizers.l2(reg)
    )(s_layer_4)

    return keras.Model(inputs=input, outputs=[s_layer_5, t_layer_5])

# Defining the type of normalizing flows model used. This model uses real-valued non-volume preserving (RealNVP) transformations
class RealNVP(keras.Model):
    """Normalizing flow emulator based on Real-valued Non-Volume Preserving (RealNVP) transformations.

    Implements a RealNVP normalizing flow (Dinh et al. 2017) as a Keras model.
    The flow is composed of *num_coupling_layers* affine coupling layers whose
    parameters are predicted by the :func:`Coupling` networks.  Alternating
    binary masks ensure that every input dimension is transformed by at least
    one coupling layer.

    The model is trained to maximize the log-likelihood of the training data
    under a six-dimensional isotropic Gaussian base distribution.

    Parameters
    ----------
    num_coupling_layers : int
        Number of coupling layers in the flow.  Must be even so that the
        alternating masks cover all dimensions equally.
    """

    def __init__(self, num_coupling_layers):
        """Initialize the RealNVP model.

        Parameters
        ----------
        num_coupling_layers : int
            Number of coupling layers to stack.
        """
        super(RealNVP, self).__init__()

        self.num_coupling_layers = num_coupling_layers

        self.distribution = tfp.distributions.MultivariateNormalDiag(
            loc=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], scale_diag= [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        self.masks = np.array(
            [[1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 1], ] * (num_coupling_layers // 2), dtype="float32"
        )
        self.loss_tracker = keras.metrics.Mean(name="loss")
        self.layers_list = [Coupling(6) for i in range(num_coupling_layers)]

    @property
    def metrics(self):
        """List of the model's metrics.
        We make sure the loss tracker is listed as part of `model.metrics`
        so that `fit()` and `evaluate()` are able to `reset()` the loss tracker
        at the start of each epoch and at the start of an `evaluate()` call.
        """
        return [self.loss_tracker]

    def call(self, x, training=True):
        """Apply the normalizing flow transformation.

        Passes *x* through all coupling layers in the appropriate direction:
        forward (data → latent) when ``training=True``, or inverse (latent →
        data) when ``training=False``.

        Parameters
        ----------
        x : tf.Tensor
            Input tensor of shape ``(batch, 6)``.
        training : bool, optional
            If ``True`` (default), apply the forward (encoding) transformation.
            If ``False``, apply the inverse (decoding) transformation.

        Returns
        -------
        x : tf.Tensor
            Transformed tensor of shape ``(batch, 6)``.
        log_det_inv : tf.Tensor
            Accumulated log absolute determinant of the Jacobian, shape
            ``(batch,)``.  Positive values indicate volume expansion.
        """
        log_det_inv = 0
        direction = 1
        if training:
            direction = -1
        for i in range(self.num_coupling_layers)[::direction]:
            x_masked = x * self.masks[i]
            reversed_mask = 1 - self.masks[i]
            s, t = self.layers_list[i](x_masked)
            s *= reversed_mask
            t *= reversed_mask
            gate = (direction - 1) / 2
            x = (
                reversed_mask
                * (x * tf.exp(direction * s) + direction * t * tf.exp(gate * s))
                + x_masked
            )
            log_det_inv += gate * tf.reduce_sum(s, [1])

        return x, log_det_inv

    def log_loss(self, data):
        """Compute the negative mean log-likelihood of a batch.

        Parameters
        ----------
        data : tf.Tensor
            Batch tensor of shape ``(batch, n_features + 1)``.  The last
            column is ignored; the first ``n_features`` columns are the
            input features.

        Returns
        -------
        tf.Tensor
            Scalar negative mean log-likelihood for the batch.
        """
        x = data[:,0:-1]
        m = data[:,0]
        y, logdet = self(x)
        log_likelihood = self.distribution.log_prob(y) + logdet
        return -tf.reduce_mean(log_likelihood)

    def train_step(self, data):
        """Perform a single gradient-descent training step.

        Parameters
        ----------
        data : tf.Tensor
            Training batch passed by Keras during ``model.fit()``.

        Returns
        -------
        dict
            Dictionary ``{"loss": <mean loss for this batch>}``.
        """
        with tf.GradientTape() as tape:
            loss = self.log_loss(data)

        g = tape.gradient(loss, self.trainable_variables)
        self.optimizer.apply_gradients(zip(g, self.trainable_variables))
        self.loss_tracker.update_state(loss)

        return {"loss": self.loss_tracker.result()}

    def test_step(self, data):
        """Perform a single evaluation step.

        Parameters
        ----------
        data : tf.Tensor
            Evaluation batch passed by Keras during ``model.evaluate()``.

        Returns
        -------
        dict
            Dictionary ``{"loss": <mean loss for this batch>}``.
        """
        loss = self.log_loss(data)
        self.loss_tracker.update_state(loss)

        return {"loss": self.loss_tracker.result()}
