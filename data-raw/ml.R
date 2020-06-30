library(magrittr)
library(keras)

work_dir <- '/home/crumplecup/work'
png_dir <- '/home/crumplecup/work/png3band'
jpg_dir <- '/home/crumplecup/work/jpg3band'
in_dir <- '/home/crumplecup/work/rip_3band'
# dir.create(png_dir)
# dir.create(jpg_dir)


files <- list.files(in_dir)
files <- files[!(1:length(files)) %in% grep('aux', files)]

setwd(png_dir)
for (file in files) {
  img <- raster::stack(file.path(in_dir, file))
  fname <- strsplit(file, '.img') %>% unlist
  png(fname)
  plotRGB(img)
  dev.off()
}

setwd(jpg_dir)
for (file in files) {
  img <- raster::stack(file.path(in_dir, file))
  fname <- strsplit(file, '.img') %>% unlist
  jpeg(paste0(fname,'.jpeg'))
  raster::plotRGB(img)
  dev.off()
}

in_dir <- '/home/crumplecup/work/rip_4band'
files <- list.files(in_dir)
files <- files[!(1:length(files)) %in% grep('aux', files)]
ls <- list()

for (i in seq_along(files)) {
  img <- raster::stack(file.path(in_dir, files[i]))
  fname <- strsplit(files[i], '.img') %>% unlist
  ar <- raster::values(img)
  ar <- ar[!is.na(ar[ ,1]), ]
  ls[[i]] <- ar
}

car <- array(0, c(length(ls), 15))
for (i in 1:length(ls)) {
  mat <- ls[[i]]
  car[i, 1] <- mean(mat[ , 1])
  car[i, 2] <- sd(mat[ , 1])
  car[i, 3] <- sum(mat[ , 1])
  car[i, 4] <- mean(mat[ , 2])
  car[i, 5] <- sd(mat[ , 2])
  car[i, 6] <- sum(mat[ , 2])
  car[i, 7] <- mean(mat[ , 3])
  car[i, 8] <- sd(mat[ , 3])
  car[i, 9] <- sum(mat[ , 3])
  car[i, 10] <- mean(mat[ , 4])
  car[i, 11] <- sd(mat[ , 4])
  car[i, 12] <- sum(mat[ , 4])
  ndvi <- (mat[ , 4] - mat[ , 1]) / (mat[ , 4] + mat[ , 1])
  car[i, 13] <- mean(ndvi)
  car[i, 14] <- sd(ndvi)
  car[i, 15] <- sum(ndvi)

}

# create labels from filenames
labs <- vector(length(files), mode = 'numeric')
labs[grep('pc', files)] <- 1
labs[grep('fc', files)] <- 2

df <- data.frame(cov = labs,
                 red_mn = car[ , 1],
                 red_sd = car[ , 2],
                 red_sm = car[ , 3],
                 grn_mn = car[ , 4],
                 grn_sd = car[ , 5],
                 grn_sm = car[ , 6],
                 blu_mn = car[ , 7],
                 blu_sd = car[ , 8],
                 blu_sm = car[ , 9],
                 nir_mn = car[ , 10],
                 nir_sd = car[ , 11],
                 nir_sm = car[ , 12],
                 ndv_mn = car[ , 13],
                 ndv_sd = car[ , 14],
                 ndv_sm = car[ , 15]
                 )

setwd(work_dir)
save(df, file = 'rip_4band_df.rds')

df$cov1 <- df$cov
df$cov1[df$cov == 2] <- 1

dt <- df
dt <- dt[ , c(17, 2:16)]


begin <- Sys.time()
mods <- stat_tab(dt, type = 'binom')
end <- Sys.time()
end - begin

bin1 <- mods[rev(order(mods$R2)),]
rip_bin1 <- lm(as.character(bin1[1,1]), data = dt)

dt <- df
dt <- dt[ , -17]
dt <- dt[dt$cov > 0, ]
dt$cov <- dt$cov - 1

begin <- Sys.time()
mods <- stat_tab(dt, type = 'binom')
end <- Sys.time()
end - begin

bin2 <- mods[rev(order(mods$R2)),]
rip_bin2 <- lm(as.character(bin2[1,1]), data = dt)

dt <- df
dt <- dt[ , -17]

begin <- Sys.time()
mods <- stat_tab(dt)
end <- Sys.time()
end - begin

lmod <- mods[rev(order(mods$R2)),]
save(bin1, file = 'rip_bin1.rds')
save(bin2, file = 'rip_bin2.rds')
save(lmod, file = 'rip_lmod.rds')

labs <- vector(length(files), mode = 'numeric')
labs[grep('pc', files)] <- 1
labs[grep('fc', files)] <- 2

load('rip_4band_df.rds')
load('rip_bin1.rds')
load('rip_bin2.rds')
rip_lmod <- lm(as.character(lmod[1,1]), data = df)
rip_lmod3 <- lm(as.character(lmod[1,1]), data = df)

setwd('/home/crumplecup/work/muddier')
usethis::use_data(rip_lmod)
usethis::use_data(rip_bin1)
usethis::use_data(rip_bin2)
usethis::use_data(rip_lmod3)


img_to_array <- function(files, in_path)  {
  rls <- 1:length(files) %>% as.list
  nr <- 0
  nc <- 0

  for (i in 1:length(files)) {
    rls[[i]] <- raster::stack(file.path(in_path, files[i]))
    nr[i] <- nrow(rls[[i]])
    nc[i] <- ncol(rls[[i]])
  }

  rar <- array(255, c(length(files), max(nr), max(nc), 4))
  for (i in 1:length(files)) {
    for (j in 1:4) {
      ras <- raster::raster(rls[[i]], layer = j)
      mat <- matrix(raster::values(ras), ncol = ncol(ras), byrow = TRUE)
      rar[i, 1:nrow(mat), 1:ncol(mat), j] <- mat
    }
  }
  rar[is.na(rar)] <- 255
  rar
}

begin <- Sys.time()
rip_4band <- img_to_array(files, in_dir)
end <- Sys.time()
end - begin

save(rip_4band, file = 'rip_4band.rds')
load('rip_4band.rds')
rip_4band <- list(rip_4band, labs)
save(rip_4band, file = 'rip_4band.rds')

base_dir <- '/home/crumplecup/work/ml'
dir.create(base_dir)

train_dir <- file.path(base_dir, 'train')
validation_dir <- file.path(base_dir, 'validation')
test_dir <- file.path(base_dir, 'test')
dir.create(train_dir)
dir.create(validation_dir)
dir.create(test_dir)

train_bg_dir <- file.path(train_dir, 'bg')
train_pc_dir <- file.path(train_dir, 'pc')
train_fc_dir <- file.path(train_dir, 'fc')
dir.create(train_bg_dir)
dir.create(train_pc_dir)
dir.create(train_fc_dir)

validation_bg_dir <- file.path(validation_dir, 'bg')
validation_pc_dir <- file.path(validation_dir, 'pc')
validation_fc_dir <- file.path(validation_dir, 'fc')
dir.create(validation_bg_dir)
dir.create(validation_pc_dir)
dir.create(validation_fc_dir)

test_bg_dir <- file.path(test_dir, 'bg')
test_pc_dir <- file.path(test_dir, 'pc')
test_fc_dir <- file.path(test_dir, 'fc')
dir.create(test_bg_dir)
dir.create(test_pc_dir)
dir.create(test_fc_dir)

in_dir <- png_dir
files <- list.files(in_dir) %>% sample
files_n <- length(files)
qtr_n <- ceiling(files_n / 4)
test_files <- files[1:qtr_n]
val_files <- files[(qtr_n + 1):(2*qtr_n)]
train_files <- files[(2*qtr_n + 1):files_n]

fnames <- train_files[grep('bg', train_files)]
file.copy(file.path(in_dir, fnames), file.path(train_bg_dir))
fnames <- train_files[grep('pc', train_files)]
file.copy(file.path(in_dir, fnames), file.path(train_pc_dir))
fnames <- train_files[grep('fc', train_files)]
file.copy(file.path(in_dir, fnames), file.path(train_fc_dir))

fnames <- val_files[grep('bg', val_files)]
file.copy(file.path(in_dir, fnames), file.path(validation_bg_dir))
fnames <- val_files[grep('pc', val_files)]
file.copy(file.path(in_dir, fnames), file.path(validation_pc_dir))
fnames <- val_files[grep('fc', val_files)]
file.copy(file.path(in_dir, fnames), file.path(validation_fc_dir))

fnames <- test_files[grep('bg', test_files)]
file.copy(file.path(in_dir, fnames), file.path(test_bg_dir))
fnames <- test_files[grep('pc', test_files)]
file.copy(file.path(in_dir, fnames), file.path(test_pc_dir))
fnames <- test_files[grep('fc', test_files)]
file.copy(file.path(in_dir, fnames), file.path(test_fc_dir))






datagen <- image_data_generator(rescale = 1/255)

train_generator <- flow_images_from_directory(
  '/home/crumplecup/work/ml/train',
  datagen,
  target_size = c(150, 150),
  batch_size = 20,
  class_mode = 'categorical'
)

validation_generator <- flow_images_from_directory(
  validation_dir,
  datagen,
  target_size = c(150, 150),
  batch_size = 20,
  class_mode = 'categorical'
)

batch <- generator_next(genny)

gen_ar <- function(in_path) {
  files <- list.files(in_path)
  ar <- array(0, c(length(files), 150, 150, 3))
  labs <- vector(length(files), mode = 'numeric')
  for(i in seq_along(files)) {
    img <- image_load(file.path(in_path, files[i]), target_size = c(150, 150))
    img_ar <- image_to_array(img)
    img_ar <- array_reshape(img_ar, c(1, 150, 150, 3))

    if(length(grep('bg', file) == 1)) lab <- 0
    if(length(grep('pc', file) == 1)) lab <- 1
    if(length(grep('fc', file) == 1)) lab <- 2

    ar[i, , , ] <- img_ar
    labs[i] <- lab
  }
  return(list(ar, labs))
}

dat <- gen_ar(in_dir)

mix <- sample(1:length(files))
dat_ar <- dat[[1]]
dat_ar <- dat_ar[mix, , , ]
dat_labs <- dat[[2]]
dat_labs <- dat_labs[mix]

rip_3band <- list(dat_ar, dat_labs)
save(rip_3band, file = 'rip_3band.rds')
rm(dat, dat_ar, dat_labs)
gc()

# create labels from filenames
files <- list.files(png_dir)
labs <- vector(length(files), mode = 'numeric')
labs[grep('pc', files)] <- 1
labs[grep('fc', files)] <- 2

setwd(work_dir)
load('rip_3band.rds')
rip_3band[[2]] <- labs
save(rip_3band, file = 'rip_3band.rds')

# add ndvi layer to rip_4band

ndvi <- (dat_ar[,,,4] - dat_ar[,,,1]) / (dat_ar[,,,4] + dat_ar[,,,1])

rip_ar <- array(0, c(dim(dat_ar)[-4], 5))
rip_ar[,,,1:4] <- dat_ar
rip_ar[,,,5] <- ndvi

rip_5band <- list(rip_ar, dat_labs)
save(rip_5band, file = 'rip_5band.rds')
rm(rip_4band, dat_ar, ndvi, rip_ar)
gc()

model %>% save_model_hdf5('rip_mod4_aug1000.h5')

model %>% save_model_hdf5('rip_mod4_aug200.h5')
model %>% save_model_hdf5('rip_mod4_aug40.h5')
model %>% save_model_hdf5('rip_mod_aug40.h5')
model %>% save_model_hdf5('rip_mod_aug10.h5')

rm(train_gen, valid_gen)
gc()



model <- load_model_hdf5('/home/crumplecup/work/rip_mod4_aug200.h5')
model

img_path <- file.path('/home/crumplecup/work/png3band/')
img_nm <- list.files(img_path)
img <- image_load(file.path(img_path, img_nm[38]), target_size = c(150, 150))

img_ar <- image_to_array(img)
img_ar <- array_reshape(img_ar, c(1, 150, 150, 3))
img_ar <- img_ar / 255
plot(as.raster(img_ar[1,,,]))

layer_outputs <- lapply(model$layers[1:8], function(layer) layer$output)
activation_model <- keras_model(inputs = model$input, outputs = layer_outputs)
activations <- activation_model %>% predict(img_ar)

first_layer_activation <- activations[[1]]
dim(first_layer_activation)

plot_channel <- function(channel) {
  rotate <- function(x) t(apply(x, 2, rev))
  image(rotate(channel), axes = FALSE, asp = 1,
        col = terrain.colors(12))
}

plot_channel(first_layer_activation[1,,,4])

image_size <- 58
images_per_row <- 16

for (i in 1:8) {

  layer_activation <- activations[[i]]
  layer_name <- model$layers[[i]]$name

  n_features <- dim(layer_activation)[[4]]
  n_cols <- n_features %/% images_per_row

  png(paste0("rip_activations_", i, "_", layer_name, ".png"),
      width = image_size * images_per_row,
      height = image_size * n_cols)
  op <- par(mfrow = c(n_cols, images_per_row), mai = rep_len(0.02, 4))

  for (col in 0:(n_cols-1)) {
    for (row in 0:(images_per_row-1)) {
      channel_image <- layer_activation[1,,,(col*images_per_row) + row + 1]
      plot_channel(channel_image)
    }
  }

  par(op)
  dev.off()
}


library(keras)
model <- application_vgg16(
  weights = "imagenet",
  include_top = FALSE
)
layer_name <- "block3_conv1"
layer_name <- 'conv2d_17'
filter_index <- 1
layer_output <- get_layer(model, layer_name)$output
loss <- k_mean(layer_output[,,,filter_index])

grads <- k_gradients(loss, model$input)[[1]]

grads <- grads / (k_sqrt(k_mean(k_square(grads))) + 1e-5)
iterate <- k_function(list(model$input), list(loss, grads))
c(loss_value, grads_value) %<-%
  iterate(list(array(0, dim = c(1, 150, 150, 3))))

input_img_data <-
array(runif(150 * 150 * 3), dim = c(1, 150, 150, 3)) * 20 + 128
step <- 1
for (i in 1:40) {
  c(loss_value, grads_value) %<-% iterate(list(input_img_data))
  input_img_data <- input_img_data + (grads_value * step)
}

deprocess_image <- function(x) {
  dms <- dim(x)
  x <- x - mean(x)
  x <- x / (sd(x) + 1e-5)
  x <- x * 0.1
  x <- x + 0.5
  x <- pmax(0, pmin(x, 1))
  array(x, dim = dms)
}

generate_pattern <- function(layer_name, filter_index, size = 150) {
  layer_output <- model$get_layer(layer_name)$output
  loss <- k_mean(layer_output[,,,filter_index])
  grads <- k_gradients(loss, model$input)[[1]]
  grads <- grads / (k_sqrt(k_mean(k_square(grads))) + 1e-5)
  iterate <- k_function(list(model$input), list(loss, grads))
  input_img_data <-
  array(runif(size * size * 3), dim = c(1, size, size, 3)) * 20 + 128
  step <- 1
  for (i in 1:40) {
    c(loss_value, grads_value) %<-% iterate(list(input_img_data))
    input_img_data <- input_img_data + (grads_value * step)
  }
  img <- input_img_data[1,,,]
  deprocess_image(img)
}

library(grid)
grid.raster(generate_pattern("block3_conv1", 1))

library(grid)
library(gridExtra)
dir.create("vgg_filters")
for (layer_name in c("block1_conv1", "block2_conv1",
                     "block3_conv1", "block4_conv1")) {
  size <- 140

  png(paste0("vgg_filters/", layer_name, ".png"),
      width = 8 * size, height = 8 * size)

  grobs <- list()
  for (i in 0:7) {
    for (j in 0:7) {
      pattern <- generate_pattern(layer_name, i + (j*8) + 1, size = size)
      grob <- rasterGrob(pattern,
                         width = unit(0.9, "npc"),
                         height = unit(0.9, "npc"))
      grobs[[length(grobs)+1]] <- grob
    }
  }

  grid.arrange(grobs = grobs, ncol = 8)
  dev.off()
}

model <- application_vgg16(weights = "imagenet")
img_path <- "~/Downloads/creative_commons_elephant.jpg"
img <- image_load(img_path, target_size = c(224, 224)) %>%
image_to_array() %>%
array_reshape(dim = c(1, 224, 224, 3)) %>%
imagenet_preprocess_input()


aug_gen <- flow_images_from_data(img_ar,
                                 generator = datagen,
                                 batch_size = 1)
batch <- generator_next(aug_gen)
plot(as.raster(batch[1,,,]))

library(raster)
ras <- stack(file.path(train_bg_dir, img_path))
plotRGB(ras)



library(tfdatasets)

data_dir <- get_file(
  origin = "https://storage.googleapis.com/download.tensorflow.org/example_images/flower_photos.tgz",
  fname = "flower_photos.tgz",
  extract = TRUE
)

data_dir <- file.path(dirname(data_dir), "flower_photos")
list_ds <- file_list_dataset(file_pattern = paste0(data_dir, "/*/*"))
list_ds <- file_list_dataset(file_pattern = jpg_dir)
list_ds %>% reticulate::as_iterator() %>% reticulate::iter_next()

images <- list.files(data_dir, pattern = ".jpg", recursive = TRUE)
length(images)
classes <- list.dirs(train_dir, full.names = FALSE, recursive = FALSE)
classes

set.seed(10101)

get_label <- function(file_name, cls) {
  parts <- tf$strings$split(file_name, '_')
  parts[length(parts)] %>%
    tf$equal(cls) %>%
    tf$cast(dtype = tf$float32)
}

decode_img <- function(file_path, file_name, height = 224, width = 224) {

  size <- as.integer(c(height, width))

  file_path <- paste(file_path, file_name, sep = '/')
  file_path %>%
    tf$io$read_file() %>%
    tf$image$decode_jpeg(channels = 3) %>%
    tf$image$convert_image_dtype(dtype = tf$float32) %>%
    tf$image$resize(size = size)
}

files <- jpg_dir %>% paste0('/') %>% list.files
preprocess_path <- function(file_path, fls = files, cls = classes) {
  nm <- sample(fls, 1)
  list(
    decode_img(file_path, nm),
    get_label(nm, cls)
  )
}

labeled_ds <- list_ds %>%
  dataset_map(preprocess_path(files, classes))

labeled_ds %>%
  reticulate::as_iterator() %>%
  reticulate::iter_next()

files <- jpg_dir %>% paste0('/') %>% list.files
ds <- file_list_dataset(jpg_dir)

get_data <- function(file_dir) {
  ds <- file_list_dataset(file_dir) %>%
    dataset_map(function(record) {
    size <- as.integer(c(244L, 244L))

    record %>%
      tf$io$read_file() %>%
      tf$image$decode_jpeg(channels = 3) %>%
      tf$image$convert_image_dtype(dtype = tf$float32) %>%
      tf$image$resize(size = size) %>%
      tf$reshape(c(1L, size, 3L))


    parts <- tf$strings$split(record, '_')
    parts[length(parts),-2] %>%
      tf$equal(classes) %>%
      tf$cast(dtype = tf$float32)
    list(imgs, labs)
  })
}

ds <- get_data(jpg_dir)

batch_size <- 128
ds %>% dataset_shuffle_and_repeat(buffer_size = files_n) %>%
  dataset_batch(batch_size, drop_remainder = TRUE) %>%
  dataset_prefetch(1)

model <- keras_model_sequential() %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dense(units = 5, activation = "softmax")

model %>%
  compile(
    loss = "categorical_crossentropy",
    optimizer = "adam",
    metrics = "accuracy"
  )

model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(244, 244, 3)) %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>%
  layer_max_pooling_2d(pool_size = c(2, 2)) %>%
  layer_flatten() %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 512, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")

model %>% compile(
  loss = "categorical_crossentropy",
  optimizer = 'rmsprop',
  metrics = c("acc")
)

in_dir <- '/home/crumplecup/work/jpg3band'
files <- list.files(in_dir) %>% sample
files_n <- length(files)
qtr_n <- ceiling(files_n / 4)
test_files <- files[1:qtr_n]
val_files <- files[(qtr_n + 1):(2*qtr_n)]
train_files <- files[(2*qtr_n + 1):files_n]

train_dir <- file.path(jpg_dir, 'train')
validation_dir <- file.path(jpg_dir, 'validation')
test_dir <- file.path(jpg_dir, 'test')
dir.create(train_dir)
dir.create(validation_dir)
dir.create(test_dir)

file.copy(file.path(in_dir, train_files), file.path(train_dir))
file.copy(file.path(in_dir, val_files), file.path(validation_dir))
file.copy(file.path(in_dir, test_files), file.path(test_dir))

steps_per_epoch <- 100

history <- model %>% fit(
  get_data(train_dir),
  steps_per_epoch = steps_per_epoch,
  epochs = 10,
  validation_data = get_data(validation_dir),
  validation_steps = steps_per_epoch
)

model <- keras_model_sequential() %>%
  layer_flatten() %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dense(units = 3, activation = "softmax")

model %>%
  compile(
    loss = "categorical_crossentropy",
    optimizer = "adam",
    metrics = "accuracy"
  )


prepare <- function(ds, batch_size, shuffle_buffer_size) {

  if (shuffle_buffer_size > 0)
    ds <- ds %>% dataset_shuffle(shuffle_buffer_size)

  ds %>%
    dataset_batch(batch_size) %>%
    # `prefetch` lets the dataset fetch batches in the background while the model
    # is training.
    dataset_prefetch(buffer_size = tf$data$experimental$AUTOTUNE)
}

model %>%
  fit(
    prepare(get_data(train_dir), batch_size = 32, shuffle_buffer_size = 1000),
    epochs = 5,
    verbose = 2
  )

get_data(train_dir) %>%
  reticulate::as_iterator() %>%
  reticulate::iter_next()


conv_base <- application_vgg16(
  weights = 'imagenet',
  include_top = FALSE,
  input_shape = c(70, 71, 3)
)
conv_base

feats <- conv_base %>% predict(dat_ar[,,,1:3])
flat <- array_reshape(feats, dim = c(nrow(feats),
                                     Reduce('*', dim(feats)[-1])))

input_tensor <- layer_input(shape = dim(feats)[-1])
output_tensor <- input_tensor %>%
  layer_flatten() %>%
  layer_dense(units = 256, activation = 'relu') %>%
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 3, activation = "softmax")

model <- keras_model(input_tensor, output_tensor)

model %>%
  compile(
    loss = "categorical_crossentropy",
    optimizer = "adam",
    metrics = "accuracy"
  )

model %>% fit(feats[c(train_ids, val_ids),,,], lab_ar[c(train_ids, val_ids),],
              epochs = 10, batch_size = 128)
model %>% evaluate(feats[test_ids,,,], lab_ar[test_ids,])
model %>% save_model_hdf5('vvg_mod3_aug25.h5')



setwd('/home/crumplecup/work')
conv_mod <- load_model_hdf5('vvg_mod3_aug25.h5')
setwd('/home/crumplecup/work/muddier')
usethis::use_data(conv_mod, overwrite = T)

ras_dir <- '/media/crumplecup/BentonCo/Statewide2018_Prelim/JustBenton'
mod_path <- '/home/crumplecup/work/rip_mod5_aug50.h5'
out_dir <- ('/home/crumplecup/work')

data(samples, package = 'riparian')
data(prc_strms, package = 'riparian')
data(prc_mtls, package = 'riparian')
data(prc_per, package = 'riparian')
data(per_buff, package = 'riparian')

sam <- riparian::sample_streams(20)
sp::plot(sam[12,])

sm <- sam[c(1,3,5:12), ]

profvis::profvis(obs <- plot_samples(ras_dir, out_dir, sm[1:2, ], mod_path))
obs

setwd(out_dir)
load('sam.rds')
slc_path <- file.path(out_dir, 'slices')
thmb_path <- file.path(out_dir, 'thumbnails')
dir.create(slc_path)
dir.create(thmb_path)

# from dir of rasters, make thumbnails of sample areas and save to drive
Rprof()
thumbnails(ras_dir, thmb_path, polys = sam)
Rprof(NULL)

report <- summaryRprof()
report$by.self

data(samples, package = 'riparian')
obs <- plot_samples(ras_dir, slc_path, samples[1:3, ], mod_path)
obs

files <- get_ras(thmb_path)
model <- keras::load_model_hdf5(mod_path)
r <- raster::stack(file.path(thmb_path, files[1]))
crs_ref <- raster::crs(r)
polys <- sam
polys <- sp::spTransform(polys, crs_ref)
raster::plotRGB(r)
area <- lapply(methods::slot(polys[1,], 'polygons'),
               function(x) lapply(methods::slot(x, 'Polygons'),
                                  function(y) methods::slot(y, 'coords')))
area <- area[[1]]

for (j in 1:50)  {
  box <- riparian::spatialize(area[1], crs_ref)
  frm <- raster::extent(box)
  slc <- raster::mask(r, box)
  slc <- raster::crop(slc, frm)

  ar <- array(255, c(1, 70, 71, 5))
  for (k in 1:4) {
    ras <- raster::raster(slc, layer = k)
    mat <- matrix(raster::values(ras), ncol = ncol(ras), byrow = TRUE)
    ar[1, 1:nrow(mat), 1:ncol(mat), k] <- mat
  }
  ar[is.na(ar)] <- 255
  ar[1,,,5] <- (ar[1,,,1] - ar[1,,,4]) / (ar[1,,,1] + ar[1,,,4])
  ar[1,,,] <- ar[1,,,] / 255

  pred <- predict(model, ar)

  car <- color_array(r)
  prd <- predict(rip_lmod, newdata = car)

  out_dir <- ('/home/crumplecup/work')
  setwd(out_dir)
  load('sam.rds')
  slc_path <- file.path(out_dir, 'slices')
  ras_dir <- '/media/crumplecup/Seagate Backup Plus Drive/gis/ortho2018'
  mod_path <- '/home/crumplecup/work/mod5_simple20_2.h5'
  obs <- plot_samples(ras_dir, slc_path, sam[c(14, 1,3,5:13,16:18), ],
                      method = 'ml', mod_path)
  obs

  obs <- plot_samples(ras_dir, slc_path, sam[c(14, 1,3,5:13,16:18), ],
                      method = 'lm3')
  obs
}

