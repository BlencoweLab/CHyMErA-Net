import keras
from keras.models import Model, Input
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv1D, AveragePooling1D
from keras.optimizers import Adam
from keras.utils import multi_gpu_model


def cnn(seq_length, drop=0.22, pool_size=3, num_filters=52, filter_size=4,
        hidden_units=11, learning_rate=0.0001,
        multi_gpu=False, with_aux=True):
    seq_input = Input(shape=(seq_length, 4), name='sequence_input')
    inputs = [seq_input]

    x = Conv1D(num_filters, filter_size, activation='relu', padding='valid',
               input_shape=(seq_length, 4))(seq_input)
    x = AveragePooling1D(pool_size)(x)
    x = Dropout(drop)(x)

    x = Flatten()(x)

    if with_aux:
        aux_input = Input(shape=(with_aux,), name='aux_input')
        inputs.append(aux_input)
        x = keras.layers.concatenate([x, aux_input], axis=-1)

    x = Dense(hidden_units, activation='relu')(x)
    x = Dropout(drop)(x)

    prediction = Dense(1, activation='sigmoid')(x)

    optimizer = Adam(lr=learning_rate)
    model = Model(inputs=inputs, outputs=[prediction])

    if multi_gpu:
        model = multi_gpu_model(model, gpus=2)
    model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    return model
