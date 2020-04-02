
# keras model

setwd(work_dir)
load('rip_5band.rds')

dat_ar <- rip_5band[[1]] / 255
dat_labs <- rip_5band[[2]]
lab_ar <- to_categorical(dat_labs)
dat_n <- length(dat_labs)
qtr_n <- ceiling(dat_n / 4)
test_ids <- 1:qtr_n
val_ids <- (qtr_n + 1):(2 * qtr_n)
train_ids <- (2 * qtr_n + 1):dat_n
batch_n <- 20
train_steps <- floor(length(train_ids) / batch_n)
test_steps <- floor(length(test_ids) / batch_n)
valid_steps <- floor(length(val_ids) / batch_n)

input_tensor <- layer_input(shape = c(70, 71, 5))
output_tensor <- input_tensor %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(70, 71, 5)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 512, activation = "relu") %>%
  # layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")

model <- keras_model(input_tensor, output_tensor)
summary(model)

model %>%
  compile(
    loss = "categorical_crossentropy",
    optimizer = "adam",
    metrics = "accuracy"
  )


model %>% fit(dat_ar[train_ids,,,], lab_ar[train_ids,],
              epochs = 12, batch_size = 128)

model %>% evaluate(dat_ar[test_ids,,,], lab_ar[test_ids,])

model %>% save_model_hdf5('keras_mod5_norm10.h5')
