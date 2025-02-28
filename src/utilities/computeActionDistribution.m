function [action,action_mean,action_std] = computeActionDistribution(actor,batch_obs_tot,states_comb_norm,obs_info,n_el)
    % Action distribution computation
    n_gauss_samples = 1000;
    action = zeros(length(states_comb_norm.Variables),n_el,n_gauss_samples);
    for i = 1:length(states_comb_norm.Variables)
        batch_obs = batch_obs_tot(:,:,i);
        batch_obs = reshape(batch_obs,length(obs_info.Name),1,n_el);
        batch_obs = repmat(batch_obs,[1,1,n_gauss_samples]);
        act = actor.getAction({batch_obs});
        action(i,:,:) = reshape(act{1},n_el,n_gauss_samples);
    end
    action_mean = mean(action,3);
    action_std = std(action,0,3);
end