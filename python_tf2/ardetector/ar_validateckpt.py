# This file validates the model while training

import numpy as np
import math
import tensorflow as tf
import os
import ar_network
import ar_reader

def validate_ckpt(session,m,config,step,filemax,val_data):
    cost = 0
    TP = 0
    FP = 0
    FN = 0
    TPe = np.full(shape=(filemax), fill_value=0, dtype=np.float32)
    FPe = np.full(shape=(filemax), fill_value=0, dtype=np.float32)
    FNe = np.full(shape=(filemax), fill_value=0, dtype=np.float32)
    for batch_input, batch_target in val_data:
        #print(file)
        _,loss, TP_, FP_, FN_, TPe_, FPe_, FNe_ = predict(session,m,config,step,file)
        #print(mask)
        TP += TP_
        FP += FP_
        FN += FN_
        TPe += TPe_
        FPe += FPe_
        FNe += FNe_
        cost += sum(loss)
        #cost += tf.reduce_sum(tf.multiply(loss,mask))
        #mask_length += tf.reduce_sum(mask)
    cost = mean(cost)
    precision = TP/(TP + FP)
    recall = TP/(TP + FN)
    F1 = 2*precision*recall/(precision + recall)
    precision_e, recall_e, F1_e, threshold = CalcPerf(TPe,FPe,FNe)
    return cost, precision, recall, F1, precision_e, recall_e, F1_e, threshold

def iterate_files(config,filemax):
    filenum = 0
    val_files = open(os.path.join(config.sum_dir,'validationSet.txt'),'r')
    for file in val_files:
        #if file.endswith(".csv"):
        filenum += 1
        if filenum <= filemax:
            yield os.path.join(config.val_dir, file.replace('\n','') + '.csv')

def predict(session,m,config,step,file):
    data = ar_reader.ArousalData(file, config)
    x, t, w = next(data)
    #print(config.checkpoint_file(step))
    softmax, loss, TP_, FP_, FN_ = session.run([m.softmax, m.loss, m.TP, m.FP, m.FN], feed_dict={
        m.features: x,
   	    m.targets: t,
        m.mask: w,
        m.step: data.iter_steps,
        m.is_training: False
        })
    TPe, FPe, FNe = CalcHits(np.array([x[1] for m, x in enumerate(softmax)]),np.array([x[1] for m, x in enumerate(t)]))
    return softmax, loss, mask, TP_, FP_, FN_, TPe, FPe, FNe

def CalcPerf(TP_,FP_,FN_):
    best_precision = 0
    best_recall = 0
    best_F1 = 0
    tstep = 0.025
    precision = np.full(shape=(len(TP_)), fill_value=0, dtype=np.float32)
    recall = np.full(shape=(len(TP_)), fill_value=0, dtype=np.float32)
    F1 = np.full(shape=(len(TP_)), fill_value=0, dtype=np.float32)
    for i in range(0,len(TP_)):
        if (TP_[i] + FP_[i]) > 0:
            precision[i] = TP_[i]/(TP_[i] + FP_[i])
        else:
            precision[i] = 0
        if (TP_[i] + FN_[i]) > 0:
            recall[i] = TP_[i]/(TP_[i] + FN_[i])
        else:
            recall[i] = 0
        if (precision[i] + recall[i]) > 0:
            F1[i] = 2*precision[i]*recall[i]/(precision[i] + recall[i])
        else:
            F1[i] = 0
        if F1[i] > best_F1:
            best_F1 = F1[i]
            best_recall = recall[i]
            best_precision = precision[i]
            best_threshold = (i+1)*tstep
    return best_precision, best_recall, best_F1, best_threshold

def CalcHits(softmax,target):
    # Initialization
    tstep = 0.025
    dur = 10 + 2
    margin = 5
    thresholds = np.arange(tstep,1.0,tstep)
    TP_ = np.full(shape=(len(thresholds)), fill_value=0, dtype=np.float32)
    FP_ = np.full(shape=(len(thresholds)), fill_value=0, dtype=np.float32)
    FN_ = np.full(shape=(len(thresholds)), fill_value=0, dtype=np.float32)
    for t in range(0,len(thresholds)):
        # Threshold
        seq = 1*(softmax > thresholds[t])
        #print(target)
        target = (target == 1)
        # Connect 
        #print(seq)
        #print(target)
        for i in range(0,len(seq) - dur + 1):
            if seq[i] & seq[i + dur - 1]:
                seq[i:i+dur] = 1
        # Remove 3 sec
        ar_dur = 0
        ar = 0
        ar_start = 0
        for i in range(0,len(seq)):
            if seq[i]:
                if not ar:
                    ar = 1
                    ar_start = i
                ar_dur += 1
            elif ar & (seq[i] == 0):
                ar = 0
                if (ar_dur < 3):
                    seq[ar_start:ar_start+ar_dur] = 0
                ar_dur = 0
            if(i == (len(seq)-1)) & ar & (ar_dur < 3):
                seq[ar_start:ar_start+ar_dur] = 0
        # Calculate onset and offset'
        # Prediction
        n_ar = 0
        new_ar = True
        pred_start = np.array([])
        pred_stop = np.array([])
        for i in range(0,len(seq)):
            if seq[i] & new_ar:
                n_ar += 1
                new_ar = False
                pred_start = np.append(pred_start,i)
            if (seq[i] == 0) & (seq[i-1] == 1) & (i > 0):
                pred_stop = np.append(pred_stop,i-1)
                new_ar = True
            if (i == len(seq)-1) & seq[i]:
                pred_stop = np.append(pred_stop,i)
        # Targets
        new_ar = True
        target_start = np.array([])
        target_stop = np.array([])
        for i in range(0,len(target)):
            if target[i] & new_ar:
                new_ar = False
                target_start = np.append(target_start,i)
            if (target[i] == 0) & (target[i-1] == 1) & (i > 0):
                target_stop = np.append(target_stop,i-1)
                new_ar = True
            if (i == len(seq)-1) & target[i]:
                target_stop = np.append(target_stop,i)
        # Statistics
        target_start = target_start - margin
        target_stop = target_stop + margin
        t_scored = np.full(shape=(len(target_start)), fill_value=0, dtype=np.float32)
        p_scored = np.full(shape=(len(pred_start)), fill_value=0, dtype=np.float32)
        for i in range(0,n_ar):
            match = ((target_stop >= pred_start[i]) & (pred_stop[i] >= target_start))
            if (match).any():
                t_scored[([m for m, x in enumerate(match) if x])] = 1
                p_scored[i] = 1
            else:
                p_scored[i] = 0
        TP = sum(p_scored == 1)
        FP = sum(p_scored == 0)
        FN = sum(t_scored == 0)
        TP_[t] = TP
        FP_[t] = FP
        FN_[t] = FN
    return TP_, FP_, FN_
    

