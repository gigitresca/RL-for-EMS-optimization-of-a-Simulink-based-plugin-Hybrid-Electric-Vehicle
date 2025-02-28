function [QValue,QValueMax,QvalueMaxIdx] = computeQvalue(critic,states_comb_norm,batch_obs_tot,batch_act,n_act,n_el,obs_info)
    qvalue = zeros(length(states_comb_norm.Variables),n_act*n_el);
    for i = 1:length(states_comb_norm.Variables)
        batch_obs = batch_obs_tot(:,:,i);
        batch_obs = reshape(batch_obs,length(obs_info.Name),1,n_el);
        batch_obs = repmat(batch_obs,[1,1,n_act]);
        qvalue(i,:) = critic(1).getValue({batch_obs},{batch_act});
    end
    % Q-values reshape for plot
    QValue = zeros(n_act,n_el,length(states_comb_norm.Variables));
    QValueMax = zeros(length(states_comb_norm.Variables),n_el);
    QvalueMaxIdx = zeros(length(states_comb_norm.Variables),n_el);
    QValueFix = zeros(n_el,n_act);
    for i = 1:length(states_comb_norm.Variables)
        QValueFix = reshape(qvalue(i,:),n_el,n_act);
        [QValueMax(i,:),QvalueMaxIdx(i,:)] = max(QValueFix,[],2);
        QValue(:,:,i) = QValueFix';
    end
end