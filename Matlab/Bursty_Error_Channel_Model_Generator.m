% % Chad Cole
% % 10/11/2021
% % 
% % This is a script to create a realistic bursty error channel model
% % 
% % Assume two states:  1 Good, where PER = alpha
% %                     2 Bad,  where PER = beta
% %                     
% %                     4 State transition probabilities


function [error_out, next_state] = Bursty_Error_Channel_Model_Generator(current_state, alpha, beta, good_transition_bias)

% alpha = 0.001;
% beta = 0.1;
transition = 0.1;
% good_transition_bias = 10;

Prob_1_given_0 = transition/good_transition_bias;
Prob_0_given_1 = transition;

% % Prob_error = Prob(Good)*Prob_err|Good + Prob(Bad)*Prob_err+Bad

error_out = 0;
rand_num = rand(1);
state_rand_num = rand(1);
if (current_state == 0)
    if rand_num <= alpha
        error_out = 1;
    end
% Set up next state
    if (state_rand_num <= Prob_1_given_0)  %Transition
        next_state = 1;
    else
        next_state = current_state;
    end
else  % in state 1
    if rand_num <= beta
        error_out = 1;
    end
% Set up next state
    if (state_rand_num <= Prob_0_given_1)  %Transition
        next_state = 0;
    else
        next_state = current_state;
    end
end
