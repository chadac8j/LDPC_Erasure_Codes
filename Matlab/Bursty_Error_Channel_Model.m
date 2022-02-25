% % Chad Cole
% % 10/11/2021
% % 
% % This is a script to create a realistic bursty error channel model
% % 
% % Assume two states:  1 Good, where PER = alpha
% %                     2 Bad,  where PER = beta
% %                     
% %                     4 State transition probabilities
                    
alpha = 0.001;
beta = 0.1;
transition = 0.1;
good_transition_bias = 10;

State_transition_table = [(1 - transition/good_transition_bias) transition/good_transition_bias; transition (1 - transition)];

% % Prob_error = Prob(Good)*Prob_err|Good + Prob(Bad)*Prob_err+Bad

% % P_Good = P_Good*(1 - transition/good_transition_bias) + P_Bad*transition;
% % P_Bad + P_Good = 1;
% % P_Good/good_transition_bias = (1 - P_Good)
% % P_Good(1 + 1/good_transition_bias) = 1
% % P_Good = 1/(1 + 1/good_transition_bias);

% % Monte Carlo sim
num_trials = 1000000;
num_errors = 0;

current_state = 0;
next_state = 0;

for ii = 1:num_trials
%     cur_error = 0;
%     rand_num = rand(1);
%     state_rand_num = rand(1);
%     if (current_state == 0)
%         if rand_num <= alpha
%             cur_error = 1;
%         end
%     % Set up next state
%         if (state_rand_num <= State_transition_table(1, 2))  %Transition
%             current_state = 1;
%         end
%     else  % in state 1
%         if rand_num <= beta
%             cur_error = 1;
%         end
%     % Set up next state
%         if (state_rand_num <= State_transition_table(2, 1))  %Transition
%             current_state = 0;
%         end
%     end
    
    [cur_error, next_state] = Bursty_Error_Channel_Model_Generator(next_state, alpha, beta, good_transition_bias);
    num_errors = num_errors + cur_error;
    
end

error_rate_sim = num_errors/num_trials

% % Prob_error = Prob(Good)*Prob_err|Good + Prob(Bad)*Prob_err|Bad
alpha_vec = [0.01, 0.02, 0.03, 0.04, 0.05];
beta_vec = [0.2, 0.3, 0.4, 0.5, 0.6];

for ii = 1:length(alpha_vec)
    for jj = 1:length(beta_vec)
        disp([alpha_vec(ii), beta_vec(jj)])
        Prob_error_analytic = (1/(1 + 1/good_transition_bias))*alpha_vec(ii) + (1 - 1/(1 + 1/good_transition_bias))*beta_vec(jj)
    end
end

