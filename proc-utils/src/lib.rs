use std::str::FromStr;

use proc_macro::TokenStream;



#[proc_macro]
pub fn repeat_across_x_times(input: TokenStream) -> TokenStream {
    let input_str = input.to_string();

    let mut split = input_str.split("@");

    let expr = split.next().expect("Error parsing expresion for expand vector").trim();
    let num = split.next().expect("Error parsing number for expand vector").trim().parse::<usize>().expect("Error parsing macro number to digits");

    
    let exprs = vec![expr; num];

    let replaced_exprs = exprs
                .iter()
                .enumerate()
                .map(|(i, expr)| {
                    let index = expr.find('*').expect("Wrong expression for find across, no asterisk in it");
                    let i_string = i.to_string();

                    format! { "{}{}{}", expr[..index - 1].to_string(), i_string, expr[index + 1..].to_string() }
                })
                .collect::<Vec<String>>();
    
    let final_stream = replaced_exprs.join("\n");

    TokenStream::from_str(final_stream.as_str()).expect("Error parsing end stream")

}

#[proc_macro]
pub fn repeat_across_x_times_and_create_struct(input: TokenStream) -> TokenStream {
    let input_str = input.to_string();

    let mut split = input_str.split("@");

    let struct_pattern = split.next().expect("Error parsing struct name for expand vector").trim();
    let expr = split.next().expect("Error parsing expresion for expand vector").trim();
    let num = split.next().expect("Error parsing number for expand vector").trim().parse::<usize>().expect("Error parsing macro number to digits");

    
    let exprs = vec![expr; num];

    let replaced_exprs = exprs
                .iter()
                .enumerate()
                .map(|(i, expr)| {
                    let index = expr.find('*').expect("Wrong expression for find across, no asterisk in it");
                    let i_string = i.to_string();

                    format! { "{}{}{}", expr[..index - 1].to_string(), i_string, expr[index + 1..].to_string() }
                })
                .collect::<Vec<String>>();
    
    let prepard_stream = replaced_exprs.join("\n");
    
    let index_asterisk_pattern = struct_pattern.find('*').expect("Error with struct pattern");
    let final_stream = format!(
        "{}{}{}",
        struct_pattern[..index_asterisk_pattern].to_string(),
        prepard_stream,
        struct_pattern[index_asterisk_pattern + 1..].to_string(),
    );

    TokenStream::from_str(final_stream.as_str()).expect("Error parsing end stream")

}


#[proc_macro_attribute]
pub fn add_x_members(attr: TokenStream, item: TokenStream) -> TokenStream {
    let attr_str = attr.to_string();
    let item_str = item.to_string();

    let attr_split = attr_str.split(";");
    
    #[derive(Clone)]
    struct RepeatToken {
        pattern: String,
        number: usize,
    }

    impl RepeatToken {
        pub fn assemble(&self) -> String {
            vec![self.pattern.as_str(); self.number]
                .iter()
                .cloned()
                .enumerate()
                .map(|(i, s)| { 
                	let index = s.find('*').expect("Wrong pattern for add marker macro");
                	let i_string = i.to_string();
                	
                	format!{ "{}{}{}", s[..index - 1].to_string(), i_string, s[index + 1..].to_string() }
                 })
                .collect::<Vec<String>>()
                .join(" ")

        }
    }

    let repeat_tokens: Vec<RepeatToken> = attr_split.into_iter().enumerate().map(|(i, s)| {
        let mut s_split = s.trim().split("=>");

        let pattern = s_split
                .next()
                .expect(format!("Error parsing pattern #{i}").as_str())
                .trim()
                .to_string();
        let number: usize = s_split
                .next()
                .expect(format!("Error parsing number #{i}").as_str())
                .trim()
                .parse()
                .expect(format!("Error getting digit for {i}").as_str());

        RepeatToken { pattern , number }
    })
    .collect();   

    let joined_items = repeat_tokens
                    .iter()
                    .cloned()
                    .map(|x| x.assemble())
                    .collect::<Vec<String>>()
                    .join(" ");

    let (before_curly, after_curly) = item_str.split_once("{").expect("Error split once'ing the struct");

    let final_str = format! {
        r#"
            {before_curly}{{ {joined_items} {after_curly}
        "#
    };
	
    TokenStream::from_str(final_str.as_str()).expect("Error finalizing token stream for x members add macro")
}