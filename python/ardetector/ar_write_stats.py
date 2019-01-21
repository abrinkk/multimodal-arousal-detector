from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
import numpy as np
import os

path = 'C:\\tmp\\train_stats\\W_NORM_SW_test_c9_lr_0.0010000,kp_1.00,bs_300,nh_128,rl_20'
path_out = 'D:\\ardetector\\data\\analytics'
event_acc = EventAccumulator(path)
event_acc.Reload()

# Save train loss
path_out2 = str.join('\\',[path_out, 'train_loss.txt'])
out_file = open(path_out2,'w')
np.savetxt(out_file,event_acc.Scalars('train_loss'),delimiter=',',fmt='%.5f')
out_file.close()

# Save test loss
path_out2 = str.join('\\',[path_out, 'test_loss.txt'])
out_file = open(path_out2,'w')
np.savetxt(out_file,event_acc.Scalars('Test_cost'),delimiter=',',fmt='%.5f')
out_file.close()

# Save test f1 ar
path_out2 = str.join('\\',[path_out, 'test_f1ar.txt'])
out_file = open(path_out2,'w')
np.savetxt(out_file,event_acc.Scalars('Test_sample_F1_score'),delimiter=',',fmt='%.5f')
out_file.close()

# Save test f1 w
path_out2 = str.join('\\',[path_out, 'test_f1w.txt'])
out_file = open(path_out2,'w')
np.savetxt(out_file,event_acc.Scalars('Test_sample_F1_score_wake'),delimiter=',',fmt='%.5f')
out_file.close()

# Save train f1 ar
path_out2 = str.join('\\',[path_out, 'train_f1ar.txt'])
out_file = open(path_out2,'w')
np.savetxt(out_file,event_acc.Scalars('train_F1_ar'),delimiter=',',fmt='%.5f')
out_file.close()

# Save test f1 ar
path_out2 = str.join('\\',[path_out, 'train_f1w.txt'])
out_file = open(path_out2,'w')
np.savetxt(out_file,event_acc.Scalars('train_F1_w'),delimiter=',',fmt='%.5f')
out_file.close()