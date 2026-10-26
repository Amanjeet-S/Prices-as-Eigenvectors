function [shares, endowments, goods] = exampleEconomy()
%EXAMPLEECONOMY The three-household, three-good economy used throughout.
%   [shares, endowments, goods] = EXAMPLEECONOMY() returns the preference
%   weights and endowments of three households trading food, clothing and
%   fuel. Each household owns a different mix of goods and prefers a
%   different one, so all three markets are needed.

shares = [0.5 0.3 0.2      % household 1 leans towards food
          0.2 0.6 0.2      % household 2 leans towards clothing
          0.3 0.3 0.4];    % household 3 leans towards fuel

endowments = [4 1 0
              0 3 2
              1 0 5];

goods = ["Food"; "Clothing"; "Fuel"];
end
