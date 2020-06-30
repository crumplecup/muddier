
# keras model
library(keras)
work_dir <- '/home/crumplecup/work'
mod_path <- '/home/crumplecup/work/mod5_simple20_3.h5'
ras_dir <- '/media/crumplecup/Seagate Backup Plus Drive/ortho2018'
slc_path <- '/home/crumplecup/work/slices'
setwd(work_dir)
load(file.path(work_dir, 'sam.rds'))
load('rip_5band.rds')
logdir <- file.path(work_dir, 'model_log')
# dir.create(logdir)
# tensorboard(logdir)

dat_ar <- rip_5band[[1]] / 255
dat_labs <- rip_5band[[2]]
lab_ar <- to_categorical(dat_labs)
dat_n <- length(dat_labs)
qtr_n <- ceiling(dat_n / 4)
test_ids <- 1:qtr_n
val_ids <- (qtr_n + 1):(2 * qtr_n)
train_ids <- (2 * qtr_n + 1):dat_n
batch_n <- 128
train_steps <- floor(length(train_ids) / batch_n)
test_steps <- floor(length(test_ids) / batch_n)
valid_steps <- floor(length(val_ids) / batch_n)


callbacks_list <- list(
  callback_early_stopping(
    monitor = 'accuracy',
    patience = 7
  ),
  callback_model_checkpoint(
    filepath = mod_path,
    monitor = 'val_loss',
    save_best_only = TRUE
  ),
  callback_reduce_lr_on_plateau(
    monitor = 'loss',
    factor = 0.1,
    patience = 10
  ),
  callback_tensorboard(
    log_dir = logdir,
    histogram_freq = 1
  )
)


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

output_tensor <- input_tensor %>%
  layer_conv_2d(filters = 32, kernel_size = 1, activation = "relu",
                input_shape = c(70, 71, 5)) %>%
  # layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = "relu") %>%
  # layer_conv_2d(filters = 32, kernel_size = c(5,5), activation = "relu") %>%
  # layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dropout(rate = 0.25) %>%
  # layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 3, activation = 'softmax')


model <- keras_model(input_tensor, output_tensor)
summary(model)

model %>%
  compile(
    loss = "categorical_crossentropy",
    optimizer = 'rmsprop',
    metrics = "accuracy"
  )

begin <- Sys.time()
model %>% fit(
  dat_ar[train_ids,,,], lab_ar[train_ids,],
  epochs = 7,
  batch_size = 128,
  callbacks = callbacks_list,
  validation_data = list(dat_ar[val_ids,,,], lab_ar[val_ids, ])
)

model %>% evaluate(dat_ar[test_ids,,,], lab_ar[test_ids,])
end <- Sys.time()
end - begin

model %>% save_model_hdf5('mod5_simple20_4.h5')

mod_path <- '/home/crumplecup/work/mod5_simple20_4.h5'
# mod_path <- '/home/crumplecup/work/mod5_simple20_3.h5'
plot_samples(ras_dir, slc_path, sam[c(14, 1,3,5:13,16:18), ],
             method = 'ml', mod_path)

