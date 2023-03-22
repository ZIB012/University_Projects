import neuraluq as neuq
import neuraluq.variables as neuq_vars
from neuraluq.config import tf

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

def load_data():
    data = sio.loadmat("dataset/3_species.mat")
    t_test, u_test = data["t_test"], data["u_test"]
    t_u_train, u_train = data["t_train"], data["u_train"]
    t_f_train, f_train = data["t_f_train"], data["f_train"]
    noise = 0.1
    return t_u_train, u_train, t_f_train, f_train, noise, t_test, u_test


def ode_fn(t, u, mu):
    u1, u2, u3 = tf.split(u, 3, axis=-1)
    u1_t, u2_t, u3_t = tf.gradients(u1, t)[0], tf.gradients(u2, t)[0], tf.gradients(u3, t)[0]

    f1 = u1_t - u1 * (mu - 0.1 * u1 - 0.5 * u2 - 0.5 * u3)
    f2 = u2_t - u2 * (-mu + 0.5 * u1 - 0.3 * u3)
    f3 = u3_t - u3 * (-mu + 0.2 * u1 + 0.5 * u2)

    return tf.concat([f1, f2, f3], axis=-1)


@neuq.utils.timer
def Samplable(t_u_train, u_train, t_f_train, f_train, noise, layers):
    t_u1_train = t_u_train[:]
    u1_train, u2_train, u3_train = u_train[:, 0:1], u_train[:, 1:2], u_train[:, 2:3]
    # build processes
    process_u = neuq.process.Process(
        surrogate=neuq.surrogates.FNN(layers=layers),
        prior=neuq_vars.fnn.Samplable(layers=layers, mean=3, sigma=2),
    )
    process_r1 = neuq.process.Process(
        surrogate=neuq.surrogates.Identity(),
        prior=neuq_vars.const.Samplable(mean=2, sigma=1),
    )

    # build likelihood
    likelihood_u23 = neuq.likelihoods.Normal(
        inputs=t_u_train,
        targets=np.concatenate([u2_train, u3_train], axis=-1),
        processes=[process_u],
        out_dims=[[1, 2]],
        sigma=noise,
    )
    likelihood_u1 = neuq.likelihoods.Normal(
        inputs=t_u1_train,
        targets=u1_train,
        processes=[process_u],
        out_dims=[[0]],
        sigma=noise,
    )

    likelihood_f = neuq.likelihoods.Normal(
        inputs=t_f_train,
        targets=f_train,
        processes=[process_u, process_r1],
        pde=ode_fn,
        sigma=noise,
    )

    # build model
    model = neuq.models.Model(
        processes=[process_u, process_r1],
        likelihoods=[likelihood_u1, likelihood_u23, likelihood_f],
    )

    # assign and compile method
    # Change parameters to make the acceptance rate close to 0.6.
    method = neuq.inferences.HMC(
        num_samples=2000,
        num_burnin=5000,
        init_time_step=0.1,
        leapfrog_step=50,
        seed=777,
    )
    model.compile(method)

    # obtain posterior samples
    samples, results = model.run()
    print("Acceptance rate: %.3f \n" % (np.mean(results)))

    processes = [process_u, process_r1]
    return processes, samples, model


def plots(u_pred, t_test, u_test, t_u_train, u_train):
    u1_pred, u2_pred, u3_pred = np.split(u_pred, 3, axis=-1)
    u1_test, u2_test, u3_test = np.split(u_test, 3, axis=-1)
    u2_train, u3_train = u_train[:, 1:2], u_train[:, 2:3]  # training data
    t_u1_train = t_u_train[:]
    u1_train = u_train[:, 0:1]

    neuq.utils.hist(a_pred, name="value of $a$")
    neuq.utils.plot1d(
        t_u1_train,
        u1_train,
        t_test,
        u1_test,
        u1_pred[..., 0],
        title="inference over mode 1",
    )
    neuq.utils.plot1d(
        t_u_train,
        u2_train,
        t_test,
        u2_test,
        u2_pred[..., 0],
        title="inference over mode 2",
    )
    neuq.utils.plot1d(
        t_u_train,
        u3_train,
        t_test,
        u3_test,
        u3_pred[..., 0],
        title="inference over mode 3",
    )


if __name__ == '__main__':

    t_u_train, u_train, t_f_train, f_train, noise, t_test, u_test = load_data()
    layers = [1, 50, 50, 3]
    processes, samples, model = Samplable(t_u_train, u_train, t_f_train, f_train, noise, layers)
    u_pred, a_pred = model.predict(t_test, samples, processes, pde_fn=None)
    plots(u_pred, t_test, u_test, t_u_train, u_train)