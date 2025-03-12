function result = has_element_less_than(A, B)%存在对应位置的A<B则返回true
    % 判断A和B的长度是否相同
    if length(A) ~= length(B)
        error('数组A和B的长度必须相同');
    end
    
    % 比较A和B中对应位置的元素
    comparison = A < B;
    
    % 判断是否存在A中的元素小于B中对应位置的元素
    if any(comparison)
        result = true;
    else
        result = false;
    end
end