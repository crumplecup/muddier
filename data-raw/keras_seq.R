
load('rip_5band.rds')

dat_ar <- rip_5band[[1]]
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


datagen <- image_data_generator(rescale = 1/255)

auggen <- image_data_generator(
  rescale = 1/255,
  rotation_range = 40,
  width_shift_range = 0.2,
  height_shift_range = 0.2,
  shear_range = 0.2,
  zoom_range = 0.2,
  horizontal_flip = TRUE,
  fill_mode = 'nearest'
)


train_gen <- flow_images_from_data(
  dat_ar[train_ids,,,],
  lab_ar[train_ids,],
  auggen,
  20)

valid_gen <- flow_images_from_data(
  dat_ar[val_ids,,,],
  lab_ar[val_ids,],
  datagen,
  20)

test_gen <- flow_images_from_data(
  dat_ar[test_ids,,,],
  lab_ar[test_ids,],
  datagen,
  20)



# batch <- generator_next(valid_gen)

model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(70, 71, 5)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  # layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  # layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 512, activation = "relu") %>%
  # layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")


model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(70, 71, 5)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_dropout(rate = 0.25) %>%
  layer_conv_2d(filters = 128, kernel_size = c(5, 5), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_dropout(rate = 0.25) %>%
  layer_conv_2d(filters = 512, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_dropout(rate = 0.25) %>%
  layer_conv_2d(filters = 512, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_dropout(rate = 0.25) %>%
  layer_flatten() %>%
  layer_dense(units = 256, activation = "relu") %>%
  layer_dropout(rate = 0.25) %>%
  layer_dense(units = 512, activation = "relu") %>%
  layer_dropout(rate = 0.25) %>%
  layer_dense(units = 3, activation = "softmax")

model <- keras_model_sequential() %>%
  layer_seperable_conv_2d(filters = 32, kernel_size = 1, activation = "relu",
                input_shape = c(70, 71, 5)) %>%
  layer_seperable_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_seperable_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  # layer_dropout(rate = 0.25) %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")




summary(model)

model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = 'rmsprop',
  metrics = c("acc")
)

begin <- Sys.time()
history <- model %>% fit_generator(
  train_gen,
  steps_per_epoch = train_steps,
  epochs = 7,
  validation_data = valid_gen,
  validation_steps = valid_steps
)
model %>% evaluate_generator(test_gen, steps = test_steps)
end <- Sys.time()
end - begin


model %>% save_model_hdf5('mod5_simple20_3.h5')
model %>% save_model_hdf5('mod5_simple20_2.h5')
model %>% save_model_hdf5('mod5_simple20_1.h5')
model %>% save_model_hdf5('mod5_aug20_1.h5')
model %>% save_model_hdf5('rip_mod5_aug50.h5')

mod_path <- '/home/crumplecup/work/mod5_simple20_3.h5'
plot_samples(ras_dir, slc_path, sam[c(14, 1,3,5:13,16:18), ],
                    method = 'ml', mod_path)


