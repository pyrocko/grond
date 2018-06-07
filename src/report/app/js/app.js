'use strict';

var REPORT_APP_VERSION = 'v0.2';

function copy_properties(source, target) {
    for (var prop in source) {
        if (source.hasOwnProperty(prop)) {
            if (source[prop] != null) {
                target[prop] = source[prop];
            }
        }
    }
}

function Dummy(obj) {
    copy_properties(obj, this);
}

function ReportIndexEntry(obj) {
    copy_properties(obj, this);
    this.get_event = function() {
        if (this.event_best != null) {
            return this.event_best;
        } else {
            return this.event_reference;
        }
    }

    this.get_solution_type = function() {
        if (this.event_best != null) {
            return 'best';
        } else {
            return 'ref';
        }
    }
}

var yaml_type_map = [
    ['!grond.ReportIndexEntry', ReportIndexEntry],
    ['!pf.Event', Dummy],
    ['!pf.MomentTensor', Dummy],
    ['!grond.ParameterStats', Dummy],
    ['!grond.TargetBalancingAnalyserResult', Dummy],
    ['!grond.ResultStats', Dummy],
    ['!grond.WaveformMisfitTarget', Dummy],
    ['!grond.WaveformMisfitConfig', Dummy],
    ['!grond.WaveformTargetGroup', Dummy],
    ['!grond.GNSSCampaignMisfitTarget', Dummy],
    ['!grond.GNSSCampaignMisfitConfig', Dummy],
    ['!grond.GNSSCampaignTargetGroup', Dummy],
    ['!grond.SatelliteMisfitTarget', Dummy],
    ['!grond.SatelliteMisfitConfig', Dummy],
    ['!grond.SatelliteTargetGroup', Dummy],
    ['!grond.PhaseRatioTarget', Dummy],
    ['!grond.PhaseRatioTargetGroup', Dummy],
    ['!grond.FeatureMeasure', Dummy],
    ['!grond.CMTProblem', Dummy],
    ['!grond.RectangularProblem', Dummy],
    ['!pf.MTSource', Dummy],
    ['!pf.RectangularSource', Dummy],
    ['!pf.HalfSinusoidSTF', Dummy],
    ['!grond.PlotCollection', Dummy],
    ['!grond.PlotGroup', Dummy],
    ['!grond.PlotItem', Dummy],
    ['!grond.PNG', Dummy],
    ['!grond.PDF', Dummy],
];

function make_constructor(type) {
    var type = type;
    var construct = function(data) {
        return new type(data);
    };
    return construct;
}

var yaml_types = [];
for (var i=0; i<yaml_type_map.length; i++) {
    var type = yaml_type_map[i][1]
    var t = new jsyaml.Type(yaml_type_map[i][0], {
        kind: 'mapping',
        instanceOf: type,
        construct: make_constructor(type),
    });
    yaml_types.push(t);
}

var report_schema = jsyaml.Schema.create(yaml_types);

function parse_fields_float(fields, input, output, error, factor) {
    parse_fields(fields, input, output, error, factor, parseFloat);
}

function parse_fields_int(fields, input, output, error) {
    parse_fields(fields, input, output, error, 1.0, parseInt);
}

function parse_fields(fields, input, output, error, factor, parse) {
    for (var i=0; i<fields.length; i++) {
        var field = fields[i];
        if (input[field].length == 0) {
            val = null;
        } else {
            var val = parse(input[field]) * factor;
            if (val.isNaN) {
                error[field] = true;
                return false;
            }
        }
        output[field] = val;
    }
}


angular.module('reportApp', ['ngRoute'])

    .config(function($routeProvider, $locationProvider) {
        $locationProvider.hashPrefix('');
        $routeProvider
            .when('/', {
                controller: 'ReportController',
                templateUrl: 'report.html',
            })
            .otherwise({
                template: '<div class="container">Not found!</div>',
            });
    })

    .factory('YamlDoc', function($http) {

        var funcs = {};
        funcs.query = function(path, loaded, options) {
            $http.get(path, {'responseType': 'text'}).then(
                function(response) {
                    var doc = jsyaml.safeLoad(response.data, options);
                    loaded(doc);
                }
            );
        };
        return funcs;
    })

    .factory('YamlMultiDoc', function($http) {

        var funcs = {};
        funcs.query = function(path, loaded, options) {
            $http.get(path, {'responseType': 'text'}).then(
                function(response) {
                    jsyaml.safeLoadAll(response.data, loaded, options);
                }
            );
        };
        return funcs;
    })

    .controller('NavigationController', function($scope, $route, YamlDoc, YamlMultiDoc, $routeParams, $location) {
        $scope.$route = $route;
        $scope.$location = $location;
        $scope.$routeParams = $routeParams;

        $scope.active = function(path) {
            return (path === $location.path().substr(0,path.length)) ? 'active' : '';
        };

    })

    .factory('ReportList', function(YamlMultiDoc) {
        var report_entries = [];
        var istate = 0;
        var selected_problem_names = [];
        var funcs = {}

        funcs.reload = function() {
            report_entries.length = 0;
            selected_problem_names.length = 0;
            YamlMultiDoc.query(
                'report_list.yaml',
                function(doc) { report_entries.push(doc); istate += 1; },
                {schema: report_schema});
        };

        funcs.selected_add = function(problem_name) {
            if (!selected_problem_names.includes(problem_name)) {
                selected_problem_names.push(problem_name);
            }
        };
        
        funcs.selected_add_and_close_modal = function(problem_name) {

            funcs.selected_add(problem_name);
            $("#report-list-modal").modal('hide');
        };

        funcs.selected_remove = function(problem_name) {
            var i = selected_problem_names.indexOf(problem_name);
            if (i > -1) {
                selected_problem_names.splice(i, 1);
            }
        };

        funcs.selected_all = function() {
            selected_problem_names.length = 0;
            for (var i=0; i<report_entries.length; i++) {
                selected_problem_names.push(
                    report_entries[i].problem_name);
            }
            selected_problem_names.sort()
        };
        funcs.selected_none = function() {
            selected_problem_names.length = 0;
        };

        funcs.get_selected = function() {
            return selected_problem_names;
        };

        funcs.get_report_entries = function() {
            return report_entries;
        };

        funcs.get_istate = function() {
            return istate;
        };

        funcs.get_path = function(problem_name) {
            for (var i=0; i<report_entries.length; i++) {
                var entry = report_entries[i];
                if (entry.problem_name == problem_name) {
                    return entry.path;
                }
            }
            return null;
        };

        return funcs;
    })

    .controller('ReportListController', function($scope, YamlMultiDoc, ReportList) {
        var rl = ReportList;
        $scope.rl = rl;
        rl.reload();

        var get_order_key_funcs = {
            'event': function(x) {
                return x.get_event().name;
            },
            'time': function(x) {
                return x.get_event().time;
            },
            'magnitude': function(x) {
                return x.get_event().magnitude;
            },
            'depth': function(x) {
                return x.get_event().depth;
            },
            'type': function(x) {
                return x.get_solution_type() + '.' + x.get_event().name;
            },
            'run': function(x) {
                return x.problem_name;
            }};

        var order_skey = 'event';
        var get_order_key = get_order_key_funcs[order_skey];

        $scope.set_order = function(skey) {
            order_skey = skey;
            get_order_key = get_order_key_funcs[skey];
        };

        var ordered_lines = {};
        var ordered_lines_istate = -1;

        $scope.get_ordered_report_entries = function() {
            if (ordered_lines_istate != rl.get_istate()) {
                ordered_lines = {};
            }
            ordered_lines_istate = rl.get_istate();
            if (order_skey in ordered_lines === false) {
                var lines = [];
                for (var i=0; i<rl.get_report_entries().length; i++) {
                    var line = {
                        bg_class: null,
                        report_entry: rl.get_report_entries()[i]};

                    lines.push(line);
                }

                lines.sort(function(a, b) { 
                    var a_key = get_order_key(a.report_entry);
                    var b_key = get_order_key(b.report_entry);

                    if (a_key < b_key)
                        return -1;
                    if (a_key > b_key)
                        return 1;
                    return 0;
                });

                var i_bg_class = 0;
                var name = '';
                var bg_classes = ['even', 'odd'];
                for (var i=0; i<lines.length; i++) {
                    var rname = lines[i].report_entry.get_event().name;
                    if (rname != name) {
                        i_bg_class++;
                        name = rname;
                    }
                    lines[i].bg_class = bg_classes[i_bg_class % 2];
                }
                ordered_lines[order_skey] = lines;
            }
            return ordered_lines[order_skey];
        };

    })

    .controller('ReportController', function(
            $scope, YamlDoc, YamlMultiDoc, $timeout, $filter, $http, $sce, ReportList) {

        var groups = {};
        $scope.groups_selected = function() { return []; };
        var rl = ReportList;
        $scope.rl = rl;
        $scope.primary_problem = null;
        $scope.secondary_problem = null;
        
        $scope.compare_mode = false;

        $scope.toggle_compare_mode = function() {
            $scope.compare_mode = ! $scope.compare_mode;
            if ($scope.compare_mode) {
                for (var i=0; i<rl.get_selected().length; i++) {
                    if (rl.get_selected()[i] != $scope.primary_problem) {
                        $scope.secondary_problem = rl.get_selected()[i];
                        break;
                    }
                }
            } else {
                $scope.secondary_problem = null;
            } 
        };

        $scope.open_modal = function() {
            $('#report-list-modal').modal('show');
        };

        $scope.close_modal = function() {
            $('#report-list-modal').modal('hide');
        };

        $scope.get_selected = function() {
            var sel = rl.get_selected();
            if (sel.length == 0) {
                $scope.open_modal();
            }
            if (sel.indexOf($scope.primary_problem) == -1) {
                $scope.primary_problem = null;
            }
            if (sel.indexOf($scope.secondary_problem) == -1) {
                $scope.secondary_problem = null;
            }

            if ($scope.primary_problem === null && sel.length > 0) {
                $scope.primary_problem = sel[0];
            }
            if ($scope.secondary_problem === null && sel.length > 1) {
                $scope.secondary_problem = sel[1];
            }
            return sel;
        };

        $scope.set_primary_problem = function(problem_name) {
            $scope.primary_problem = problem_name;
        };

        $scope.set_secondary_problem = function(problem_name) {
            $scope.secondary_problem = problem_name;
        };

        var get_path = function(problem_name) {
            return rl.get_path(problem_name);
        };

        var plot_group_path = function(problem_name, group_ref) {
            return get_path(problem_name) + '/plots/' + group_ref[0] + '/' + group_ref[1] + '/' + group_ref[0] + '.' + group_ref[1];
        };

        $scope.image_path = function(group, item) {
            return plot_group_path(group.problem_name, [group.name, group.variant]) + '.' + item.name + '.d100.png';
        };

        var doc_order = function(doc1, doc2) {
            if (doc1.section < doc2.section) return -1;
            if (doc1.section > doc2.section) return 1;
            if (doc1.name < doc2.name) return -1;
            if (doc1.name > doc2.name) return 1;
            return 0;
        };

        var insert_group = function(problem_name, doc) {
            doc['problem_name'] = problem_name;
            groups[problem_name].push(doc);
            groups[problem_name].sort(doc_order);
        };

        var load_group = function(problem_name, group_ref) {
            YamlDoc.query(
                plot_group_path(problem_name, group_ref) + '.plot_group.yaml',
                function(doc) {
                    doc.template = 'group-plots';
                    insert_group(problem_name, doc);
                },
                {schema: report_schema});
        };

        var load = function(problem_name) {
            if (problem_name === null) {
                return;
            }
            if (!groups.hasOwnProperty(problem_name)) {
                console.log('load', problem_name, get_path(problem_name));
                if ($scope.primary_problem === null) {
                    $scope.primary_problem = problem_name;
                }
                groups[problem_name] = [];
                YamlDoc.query(
                    get_path(problem_name) + '/stats.yaml',
                    function(doc) {
                        doc.name = 'parameter results';
                        doc.section = 'run';
                        doc.feather_icon = 'book';
                        doc.template = 'parameter-table';

                        insert_group(problem_name, doc);
                    },
                    {schema: report_schema});

                YamlMultiDoc.query(
                    get_path(problem_name) + '/plots/plot_collection.yaml',
                    function(doc) {
                        for (var i = 0, len = doc.group_refs.length; i < len; i++)
                            load_group(problem_name, doc.group_refs[i]);
                    },
                    {schema: report_schema}
                );

                $http.get(get_path(problem_name) + '/config.yaml', {'responseType': 'text'}).then(function(data) {
                    var doc = new Dummy({
                        'name': 'config',
                        'section': 'run',
                        'feather_icon': 'code',
                        'template': 'config-file',
                        'raw_config': data.data
                    });

                    insert_group(problem_name, doc);
                });
            }
        };

        $scope.select_groups_all = function() {
            $scope.groups_selected = function() {
                return $scope.get_groups_avail();
            };
        };

        $scope.select_group_by_name = function(group_name) {
            $scope.groups_selected = function() {
                return $scope.get_groups_avail().filter(
                    function(group) {
                        if(group.name == group_name)
                            return true;
                        return false;
                    });
            };
        };

        $scope.select_group_by_section_name = function(section_name) {
            $scope.groups_selected = function() {
                return $scope.get_groups_avail().filter(
                    function(group) {
                        if(group.section == section_name)
                            return true;
                        return false;
                    });
            };
        };


        var groups_joined_cache = [];
        var cached = function(o, equals) {
            for (var i=0; i<groups_joined_cache.length; i++) {
                var x = groups_joined_cache[i];
                if (angular.equals(o, x)) {
                    return x;
                }
            }
            if (groups_joined_cache.unshift(o) > 1) {
                groups_joined_cache.length = 1;
            }
            return o;
        }

        $scope.have_groups_avail = function () {
            load($scope.primary_problem);
            return groups.hasOwnProperty($scope.primary_problem)
        };

        var tid = null;

        $scope.get_groups_avail = function() {

            load($scope.primary_problem);
            if ($scope.compare_mode) {
                load($scope.secondary_problem);
            }

            var pgroup = groups.hasOwnProperty($scope.primary_problem)
                ? groups[$scope.primary_problem] : [];

            var sgroup = groups.hasOwnProperty($scope.secondary_problem)
                ? groups[$scope.secondary_problem] : [];
            
            var group_map = new Map();
            var k = function(doc) {
                return doc.section + ' ' + doc.name;
            }
            for (var i=0; i<pgroup.length; i++) {
                group_map.set(k(pgroup[i]), [pgroup[i], null]);
            }
            if ($scope.compare_mode) {
                for (var i=0; i<sgroup.length; i++) {
                    if (group_map.has(k(sgroup[i]))) {
                        group_map.get(k(sgroup[i]))[1] = sgroup[i];
                    } else {
                        group_map.set(k(sgroup[i]), [null, sgroup[i]]);
                    }
                }
            }
            var groups2 = Array.from(group_map.values());
            var groups_joined = [];
            for (var i=0; i<groups2.length; i++) {
                var g = groups2[i][0] !== null ? groups2[i][0] : groups2[i][1];
                 groups_joined.push({
                    'section': g.section,
                    'name': g.name,
                    'template': g.template,
                    'feather_icon': g.feather_icon,
                    'can_compare': ! $scope.compare_mode || (groups2[i][0] !== null && groups2[i][1] !== null),
                    'first_in_section': true,
                    'pri': groups2[i][0],
                    'sec': groups2[i][1]});
            }

            groups_joined.sort(doc_order);

            for (var i=1; i<groups_joined.length; i++) {
                groups_joined[i].first_in_section = 
                    groups_joined[i-1].section !== groups_joined[i].section;
            }

            var c_groups_joined = cached(groups_joined);

            // don't use angular's timeout function here (infdig)!
            if (tid !== null) clearTimeout(tid);
            tid = setTimeout(function() { tid = null; feather.replace(); }, 20.)

            return c_groups_joined;
        };

        $scope.unfold_docs = function(group) {
            if ($scope.compare_mode) {
                return [group.pri, group.sec];
            } else {
                return [group.pri];
            }
        };

        $scope.group_matches_keyword = function(group) {
            if ($scope.keyword === null || $scope.keyword === undefined || $scope.keyword === '') {
                return true;
            }
            var r1;
            if (group.pri && group.pri.items) {
                r1 = $filter('filter')(group.pri.items, $scope.keyword);
            } else {
                r1 = null;
            }

            if ($scope.compare_mode) {
                var r2;
                if (group.sec && group.sec.items) {
                    r2 = $filter('filter')(group.sec.items, $scope.keyword);
                } else {
                    r2 = null;
                }

                return (r1 && r1.length > 0)
                    || (r2 && r2.length > 0);
            } else {
                return (r1 && r1.length > 0);
            }
        };

        $scope.doc_matches_keyword = function(doc) {
            return ($filter('filter')(doc.items, $scope.keyword)).length > 0;
        };
    })


    .filter('eround', function() {
        return function(input, std) {
            if (input === null || input === undefined) {
                return '-';
            } else if (std > 0) {
                var ndig = - Math.floor(Math.log10(std)) + 1;
                var factor = Math.pow(10, ndig);
                return Math.round(input * factor) / factor;
            } else {
                return input;
            }
        };
    })

    .filter('dotalign', function () {
        return function(input) {
            input = input.toString();
            var dotpos = input.indexOf('.');
            if (dotpos == -1) {
                dotpos = input.length;
            }
            var fill = ' ';
            return fill.repeat(Math.max(0, 5 - dotpos)) + input;
        };
    })

    .filter('unsafe', function($sce)
        { return $sce.trustAsHtml; })


    .filter('nounderscore', function () {
        return function (value) {
            return (!value) ? '' : value.replace(/_/g, ' ');
            };
    })

    .run(function($rootScope, $location, $anchorScroll, $routeParams) {
      //when the route is changed scroll to the proper element.
      $rootScope.$on('$routeChangeSuccess', function(newRoute, oldRoute) {
        $location.hash($routeParams.scrollTo);
        $anchorScroll();  
      });
    })

    .directive('prism', [function() {
            return {
                restrict: 'A',
                link: function ($scope, element, attrs) {
                    element.ready(function() {
                        Prism.highlightElement(element[0], Prism.languages.yaml, 'YAML');
                    });
                }
            }
    }])

    .directive('groupPlots', function() {
        return {
            templateUrl: 'templates/group_plots.tmpl.html'
        };
    })

    .directive('parameterTable', function() {
        return {
            templateUrl: 'templates/parameter_table.tmpl.html'
        };
    })

    .directive('configFile', function() {
        return {
            templateUrl: 'templates/config_file.tmpl.html'
        };
    })

    .directive('reportListModal', function() {
        return {
            templateUrl: 'templates/report_list_modal.tmpl.html'
        };
    })

    .controller('Info', function(
            $scope) {
        $scope.version = REPORT_APP_VERSION;
    });
