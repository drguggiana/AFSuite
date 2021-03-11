function cat_numbers = count_occurrences(vector_in,cat_vector)
% count the occurrences of every element of the cat_vector in the vector_in

% get the number of elements in cat_vector
num_cat = numel(cat_vector);
% allocate the output
cat_numbers = zeros(num_cat,1);
% for all the elements in the cat_vector
for el = 1:num_cat
    cat_numbers(el) = sum(vector_in==cat_vector(el));
end