### using matrix as input for training:
library(keras);
library(tensorflow);
library(ROCR)
### CNN
##change index order of arrary to the one keras package can be applied
## for 1D layers:

set.seed(123);
model = keras_model_sequential();
model %>%
  
  layer_conv_1d(filters = 75, kernel_size = 6, strides = 1,padding = "valid",input_shape = c(201,4), activation = "relu") %>%
  layer_dropout(0)%>%
  
  layer_max_pooling_1d(pool_size =4, strides =4 )%>%
  layer_dropout(0.5)%>%
  
  layer_lstm(units = 35, return_sequences = TRUE, use_bias = TRUE, activation = "relu") %>%
  layer_dropout(0.6) %>%  
  
  layer_flatten() %>%
  layer_dense(units = 2, activation = "sigmoid") %>%
  summary(model)

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = optimizer_adamax(),
  metrics = c("accuracy")
)

checkpointer = callback_model_checkpoint(filepath = "Desktop/saved_model/DeepM6A_bestmodel_201nt_R.hdf5", monitor = 'val_acc', verbose = 1, save_best_only = TRUE);
###
history <- model %>% fit(
  x_train_array_new,
  y_train,
  epochs = 70, 
  batch_size = 128,
  validation_split = 0.1,
  callbacks = checkpointer
)

model %>% evaluate(x_test_array_new, y_test, verbose = 0);

# model %>% save_model_hdf5(filepath = "Desktop/saved_model/DeepM6A_bestmodel_101nt_R.hdf5")


model = load_model_hdf5(filepath = "Desktop/saved_model/DeepM6A_bestmodel_8.hdf5")


x_inde_test = readRDS("Desktop/array_data_m6A_human/x_inde_test_100.rds");
y_inde_test = readRDS("Desktop/array_data_m6A_human/y_inde_test_100.rds");

model %>% evaluate(x_inde_test, y_inde_test, verbose = 0)

#### predict mRNA with class:
mRNAs = readRDS("Desktop/m6A_model_validataion_data/mRNA/mRNA_array_data.rds");

mRNAs_pred = model %>% predict_classes(mRNAs);
con = file("Desktop/m6A_model_validataion_data/mRNA/predicted_result.txt","w")
i = 1;
while (i<= length(mRNAs_pred))
{
  write(mRNAs_pred[i],con);
  i = i + 1
}
close(con)

#### predict with probability results:
mRNAs_pred = model %>% predict_proba(mRNAs)
mRNAs_pred_prob = as.data.frame(mRNAs_pred)
write.csv(mRNAs_pred_prob,"Desktop/m6A_model_validataion_data/mRNA/predicted_result_prob.csv",row.names = FALSE)

####lncRNAs prob:
lncRNAs=readRDS("Desktop/m6A_model_validataion_data/lncRNA/lncRNA_array_data.rds")
lncRNAs_pred = model %>% predict_proba(lncRNAs)
lncRNAs_pred_prob = as.data.frame(lncRNAs_pred)
write.csv(lncRNAs_pred_prob,"Desktop/m6A_model_validataion_data/lncRNA/predict_result_prob.csv",row.names = FALSE)

#### circRNAs prob:
circRNAs = readRDS("Desktop/m6A_model_validataion_data/circRNA/circRNA_array_data.rds")
circRNAs_pred = model %>% predict_proba(circRNAs)
circRNAs_pred_prob = as.data.frame(circRNAs_pred)
write.csv(circRNAs_pred_prob,"Desktop/m6A_model_validataion_data/circRNA/predict_result_prob.csv",row.names = FALSE)




